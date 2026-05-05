#ifndef CHIMERA_LOG_HH
#define CHIMERA_LOG_HH

#include <vector>
#include <string>
#include <utility>
#include <fstream>
#include <stdexcept>
#include <iostream>

// type for unique identifier of chimera events, can be changed to other types (e.g., std::string) if needed
using event_id_type = unsigned int;

enum class ChimeraEventType {
    GENE_FUSION_INTER_SUB_GENE_DUPLICATION,
    GENE_FUSION_INTRA_SUB_GENE_DUPLICATION,
    GENE_FUSION_EXTENDED_DELETION_FUSION,
    GENE_FUSION_EXTENDED_DELETION_REINSERTION,
    MUTATION_TRANSLOCATION,
    // both inversion breakpoints are within the same gene. this produce
    // a single chiemera record
    MUTATION_INVERSION_INTRA,
    // the two inversion breakpoints fall inside two consecutive genes
    // (in genomic order) on the same strand. this produce 2 chimera records,
    // one per affected gene: gene_A becomes head_A + RC(head_B), gene_B
    // becomes RC(tail_A) + tail_B
    MUTATION_INVERSION_INTER,
    // TODO add other types of events
};

inline const char*
event_type_to_string(const ChimeraEventType& t) {
    switch (t) {
        case ChimeraEventType::GENE_FUSION_INTER_SUB_GENE_DUPLICATION: {
            return "GENE_FUSION_INTER_SUB_GENE_DUPLICATION";
        }
        case ChimeraEventType::GENE_FUSION_INTRA_SUB_GENE_DUPLICATION: {
            return "GENE_FUSION_INTRA_SUB_GENE_DUPLICATION";
        }
        case ChimeraEventType::GENE_FUSION_EXTENDED_DELETION_FUSION: {
            return "GENE_FUSION_EXTENDED_DELETION_FUSION";
        }
        case ChimeraEventType::GENE_FUSION_EXTENDED_DELETION_REINSERTION: {
            return "GENE_FUSION_EXTENDED_DELETION_REINSERTION";
        }
        case ChimeraEventType::MUTATION_TRANSLOCATION: {
            return "MUTATION_TRANSLOCATION";
        }
        case ChimeraEventType::MUTATION_INVERSION_INTRA: {
            return "MUTATION_INVERSION_INTRA";
        }
        case ChimeraEventType::MUTATION_INVERSION_INTER: {
            return "MUTATION_INVERSION_INTER";
        }
    }
    return "UNKNOWN";
}

// a GeneLocus identifies a specific gene (or locus) in a specific genome,
// and is used to specify the source of a contribution in a chimera event
struct GeneLocus
{
    int genome_id; // id of the genome to which the gene belongs
    int gene_id; // id of the gene (or locus) in the genome

    bool operator==(const GeneLocus& other) const {
        return genome_id == other.genome_id && gene_id == other.gene_id;
    }
};

// helper function to create a GeneLocus struct from genome_id and gene_id
inline GeneLocus
make_locus(const int genome_id, const int gene_id) {
    GeneLocus locus;
    locus.genome_id = genome_id;
    locus.gene_id = gene_id;
    return locus;
}

// a chimera contribution represents a specific gene contribution to a chimera event.
// it contains all the relevant information about the contribution (e.g., source gene, offsets, length, etc.)
// it is used to keep track of the details of each contribution in a specific chimera event
struct ChimeraContribution {
    GeneLocus contribution_gene; // gene (or locus) from which the contribution comes
    int contribution_offset; // offset (int) inside the contribution gene from which the contribution starts
    int contribution_length; // length (int) of the contribution

    // snapshot of the donor gene's coordinates and strand at the moment the
    // event was applied. recorded before the mutation modifies the genome.
    int donor_start_at_event;
    int donor_end_at_event;
    int donor_strand_at_event;

    // true iff the chunk was physically reverse-complemented before being
    // inserted into the acceptor.
    //
    // set to true by:
    //   - cross-strand inter sub-gene duplication
    //   - cross-strand translocation
    //   - every contribution of an inversion event (always, even when donor
    //     and acceptor share the same strand)
    bool reverse_complemented;
#ifdef DEBUG
std::string contribution_sequence; // the actual sequence of the contribution (for debugging purposes, not used in the chimera generation process)
#endif
    std::ostream&
    operator<<(std::ostream& os) const {
        os << "gene " << contribution_gene.gene_id << " in genome " << contribution_gene.genome_id
           << " (offset " << contribution_offset << ", length " << contribution_length
           << ", donor@event start=" << donor_start_at_event
           << " end=" << donor_end_at_event
           << " strand=" << donor_strand_at_event
           << ", rc=" << (reverse_complemented ? "true" : "false") << ")";
#ifdef DEBUG
        os << ", sequence: " << contribution_sequence;
#endif
        return os;
    }
};

std::ostream&
operator<<(std::ostream& os, const ChimeraContribution& c) {
    return c.operator<<(os);
}

inline ChimeraContribution
make_chimera_contribution(
    const GeneLocus& contribution_gene, const int contribution_offset,
    const int contribution_length,
    const int donor_start_at_event, const int donor_end_at_event,
    const int donor_strand_at_event, const bool reverse_complemented
#ifdef DEBUG
    , const std::string& contribution_sequence = ""
#endif
) {
    ChimeraContribution c;
    c.contribution_gene = contribution_gene;
    c.contribution_offset = contribution_offset;
    c.contribution_length = contribution_length;
    c.donor_start_at_event = donor_start_at_event;
    c.donor_end_at_event = donor_end_at_event;
    c.donor_strand_at_event = donor_strand_at_event;
    c.reverse_complemented = reverse_complemented;
#ifdef DEBUG
    c.contribution_sequence = contribution_sequence;
#endif
    return c;
}

// a chimera acceptor represents the gene (or locus) to which the contribution is added in a chimera event
struct ChimeraAcceptor {
    GeneLocus acceptor_gene; // gene (or locus) to which the contribution is added
    // acceptor_offset is gene-relative: offset (int) inside the acceptor gene
    // measured from acceptor_gene.start.
    //
    // gene-relative offsets stay valid across subsequent mutations
    // in the same evolution step that may shift genome-absolute positions
    // (e.g. an extended-deletion "reinsertion" inserts bytes earlier in the
    // genome and would invalidate any genome-absolute offset previously
    // recorded for a gene that lies after the insertion site). The acceptor
    // gene itself moves as a block, so the relative offset within it is
    // stable.
    //
    // To recover the genome-absolute position from a chimera log entry, look
    // up the acceptor gene's post-mutation start (e.g. from the .genes file)
    // and add this offset, or use `acceptor_start_at_event` below to recover
    // the position at the time of the event
    int acceptor_offset;

    // snapshot of the acceptor gene's coordinates and strand at the moment
    // the event was applied. recorded before the mutation modifies the
    // genome (so e.g. for an inter duplication these are the target gene's
    // start/end before the inserted chunk extends `end`).
    //
    // same as the donor_*_at_event fields on ChimeraContribution:
    // subsequent mutations (in particular MUTATION_INVERSION) can re-map
    // start/end and flip strand for the same acceptor gene id; the snapshot
    // preserves the at-event state so a single record stays interpretable
    // without replaying the mutation history.
    int acceptor_start_at_event;
    int acceptor_end_at_event;
    int acceptor_strand_at_event;
};
inline ChimeraAcceptor
make_chimera_acceptor(
    const GeneLocus& acceptor_gene, const int acceptor_offset,
    const int acceptor_start_at_event, const int acceptor_end_at_event,
    const int acceptor_strand_at_event
) {
    ChimeraAcceptor a;
    a.acceptor_gene = acceptor_gene;
    a.acceptor_offset = acceptor_offset;
    a.acceptor_start_at_event = acceptor_start_at_event;
    a.acceptor_end_at_event = acceptor_end_at_event;
    a.acceptor_strand_at_event = acceptor_strand_at_event;
    return a;
}


// a chimera record represents a specific chimera event.
// it contains the acceptor gene to which the contribution is added, and a
// list of contributions that are added to the acceptor gene in the event.
struct ChimeraRecord {
    ChimeraEventType event_type; // type of the chimera event for this contribution
    ChimeraAcceptor acceptor; // the acceptor gene and its associated information to which the contribution is added
    std::vector<ChimeraContribution> contributions; // list of contributions that are added to the acceptor gene

    std::ostream&
    operator<<(std::ostream& os) const {
        os << "event type: " << event_type_to_string(event_type) << "\n";
        os << "acceptor gene: " << acceptor.acceptor_gene.gene_id << " in genome " << acceptor.acceptor_gene.genome_id
           << " (offset " << acceptor.acceptor_offset
           << ", acceptor@event start=" << acceptor.acceptor_start_at_event
           << " end=" << acceptor.acceptor_end_at_event
           << " strand=" << acceptor.acceptor_strand_at_event << ")\n";
        os << "contributions:\n";
        for (const auto& c : contributions) {
            os << "  - " << c << "\n";
        }
        return os;
    }
};

std::ostream&
operator<<(std::ostream& os, const ChimeraRecord& r) {
    return r.operator<<(os);
}

inline ChimeraRecord
make_chimera_record(
    const ChimeraEventType& event_type, const ChimeraAcceptor& acceptor,
    const std::vector<ChimeraContribution>& contributions
) {
    ChimeraRecord r;
    r.event_type = event_type;
    r.acceptor = acceptor;
    r.contributions = contributions;
    return r;
}

inline ChimeraRecord
make_chimera_record(
    const ChimeraEventType& event_type, ChimeraAcceptor&& acceptor,
    std::vector<ChimeraContribution>&& contributions
) {
    ChimeraRecord r;
    r.event_type = event_type;
    r.acceptor = std::move(acceptor);
    r.contributions = std::move(contributions);
    return r;
}



// this class keep track of all the chimeric events
class ChimeraLog {
private:

    // each entry in the log corresponds to a specific chimeric event,
    // and contains all the relevant information about the event (e.g., type of event, genome and gene involved, contribution details, etc.)
    struct ChimeraLogEntry {
        event_id_type event_id; // unique identifier for the chimera event
        ChimeraRecord record; // the record of the chimera event, containing all the details about the event
    };
    // running integer to assign unique identifiers to chimera events
    event_id_type next_event_id;

    std::vector<ChimeraLogEntry> log_entries; // vector to store all the chimera log entries


    inline void
    write_to_file(const std::string& filename, char delimiter) const {
        std::ofstream out;
        
        out.open(filename);
        if (!out) {
            throw std::runtime_error("ChimeraLog: cannot open file '" + filename + "' for writing");
        }
        std::cout << "writing " << log_entries.size() << " chimera events to '" << filename << "'\n";

        // header
        out << "event_id"                  << delimiter
            << "event_type"                << delimiter
            << "acceptor_genome_id"        << delimiter
            << "acceptor_gene_id"          << delimiter
            << "acceptor_offset"           << delimiter
            << "acceptor_start_at_event"   << delimiter
            << "acceptor_end_at_event"     << delimiter
            << "acceptor_strand_at_event"  << delimiter
            << "donor_genome_id"           << delimiter
            << "donor_gene_id"             << delimiter
            << "donor_offset"              << delimiter
            << "contribution_length"       << delimiter
            << "donor_start_at_event"      << delimiter
            << "donor_end_at_event"        << delimiter
            << "donor_strand_at_event"     << delimiter
            << "reverse_complemented"
    #ifdef DEBUG
            << delimiter << "contribution_sequence"
    #endif
            << '\n';

        // one row per contribution
        for (const auto& entry : log_entries) {
            const ChimeraRecord& r = entry.record;
            event_id_type current_event_id = entry.event_id;
            for (const auto& c : r.contributions) {
                out << current_event_id                       << delimiter
                    << event_type_to_string(r.event_type)     << delimiter
                    << r.acceptor.acceptor_gene.genome_id     << delimiter
                    << r.acceptor.acceptor_gene.gene_id       << delimiter
                    << r.acceptor.acceptor_offset             << delimiter
                    << r.acceptor.acceptor_start_at_event     << delimiter
                    << r.acceptor.acceptor_end_at_event       << delimiter
                    << r.acceptor.acceptor_strand_at_event    << delimiter
                    << c.contribution_gene.genome_id          << delimiter
                    << c.contribution_gene.gene_id            << delimiter
                    << c.contribution_offset                  << delimiter
                    << c.contribution_length                  << delimiter
                    << c.donor_start_at_event                 << delimiter
                    << c.donor_end_at_event                   << delimiter
                    << c.donor_strand_at_event                << delimiter
                    << (c.reverse_complemented ? 1 : 0)
    #ifdef DEBUG
                    << delimiter << c.contribution_sequence
    #endif
                    << '\n';
            }
        }

        if (!out) {
            throw std::runtime_error("ChimeraLog: error while writing to '" + filename + "'");
        }

        out.flush();
        out.close();
    }

public:
    ChimeraLog() : next_event_id(0), log_entries() {
        // constructor initializes the next_event_id to 0, and log_entries to an empty vector
    }
    
    inline event_id_type
    add_chimera_event(const ChimeraRecord& record) { // copy
        event_id_type id = next_event_id;
        ++next_event_id;
        log_entries.push_back(ChimeraLogEntry{id, record});
        return id;
    }
    inline event_id_type
    add_chimera_event(ChimeraRecord&& record) { // move
        event_id_type id = next_event_id;
        ++next_event_id;
        log_entries.push_back(ChimeraLogEntry{id, std::move(record)});
        return id;
    }

    inline void
    write_to_csv(const std::string& filename) const {
        write_to_file(filename, ',');
    }

    inline void
    write_to_tsv(const std::string& filename) const {
        write_to_file(filename, '\t');
    }

    inline const std::vector<ChimeraLogEntry>&
    entries() const {
        return log_entries;
    }

    inline size_t
    size() const {
        return log_entries.size();
    }
    inline bool
    empty() const {
        return log_entries.empty();
    }

};

#endif