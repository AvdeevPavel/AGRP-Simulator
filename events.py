import numpy as np


class EventType:
    all = ["rearrangement", "deletion", "insertion", "duplication"]
    REARRANGEMENT, DELETION, INSERTION, DUPLICATION = all


class RearrangementType:
    all = ["reversal", "translocation"]
    REVERSAL, TRANSLOCATION = all


class DistanceParameters:
    def __init__(self, events=0,ins_p=0, del_p=0, indel_length=0,
                 dupl_p=0, dupl_length=0, dupl_inter_p = 0,
                 num_chr=1):
        self.event_number = events
        self.rearrangement_p = 1 - del_p - ins_p - dupl_p

        self.insertion_p = ins_p
        self.deletion_p = del_p
        self.indel_length = indel_length

        self.duplication_p = dupl_p
        self.duplication_length = dupl_length
        self.duplication_interchromosome_p = dupl_inter_p
        assert np.isclose(self.rearrangement_p + self.insertion_p + self.deletion_p + self.duplication_p, 1)

        if num_chr > 1:
            self.reversal_p = .9
            self.translocation_p = .1
        else:
            self.reversal_p = 1
            self.translocation_p = 0


# Regular rearrangement operations
def apply_random_reversal(genome):
    chromosome = np.random.choice(genome.chromosomes)
    bp = sorted(np.random.choice(len(chromosome.gene_order) + 1, 2, replace=False))
    chromosome.gene_order[bp[0]:bp[1]] = reversed([-x for x in chromosome.gene_order[bp[0]:bp[1]]])


def apply_random_translocation(genome):
    chromosomes = np.random.choice(genome.chromosomes, 2, replace=False)

    bp1, bp2 = 0, 0
    while bp1 == 0 and bp2 == 0:
        bp1 = np.random.choice(len(chromosomes[0].gene_order))
        bp2 = np.random.choice(len(chromosomes[1].gene_order))

    chromosomes[0].gene_order[bp1:], chromosomes[1].gene_order[bp2:] = \
        chromosomes[1].gene_order[bp2:], chromosomes[0].gene_order[bp1:]


def apply_random_rearrangement(param, genome):
    rearrangement = np.random.choice([RearrangementType.REVERSAL, RearrangementType.TRANSLOCATION],
                                     1,
                                     p=[param.reversal_p, param.translocation_p])
    if rearrangement == RearrangementType.REVERSAL:
        apply_random_reversal(genome)
    elif rearrangement == RearrangementType.TRANSLOCATION:
        apply_random_translocation(genome)
    else:
        raise RuntimeError("Unknown rearrangement type.")


# Insertion and deletion events
def apply_random_deletion(genome, deletion_length_range):
    chromosome = np.random.choice(genome.chromosomes)
    bp = np.random.choice(chromosome.length())
    length = np.random.choice(deletion_length_range)

    if bp + length > chromosome.length():
        length = chromosome.length() - bp

    chromosome.gene_order[bp:bp + length] = []

    # remove chromosome if empty:
    if len(chromosome.gene_order) == 0:
        genome.chromosomes.remove(chromosome)


def apply_random_insertion(genome, gene, insertion_length_range):
    chromosome = np.random.choice(genome.chromosomes)
    bp = np.random.choice(chromosome.length())
    length = np.random.choice(insertion_length_range)
    block = range(gene, gene + length)
    chromosome.gene_order[bp:bp] = block
    return gene + length


# Duplication events
def apply_random_segmental_duplication(genome, duplication_length_range, interchromosome_p):
    chromosome = np.random.choice(genome.chromosomes)
    bp = np.random.choice(chromosome.length())

    length = np.random.choice(duplication_length_range)
    if bp + length > chromosome.length():
        length = chromosome.length() - bp

    position = np.random.choice(range(bp + 1) + range(bp + length, chromosome.length()))

    block = chromosome.gene_order[bp:bp + length]

    # check if duplication should be at new chromosome
    if np.random.uniform(0, 1) < interchromosome_p and len(genome.chromosomes) > 1:
        old_chr, new_chr = str(chromosome), str(chromosome)
        while old_chr == new_chr:
            chromosome = np.random.choice(genome.chromosomes)
            new_chr = str(chromosome)
        position = np.random.choice(chromosome.length())

    # apply dup:
    chromosome.gene_order[position:position] = block


def apply_wgd(genome):
    dup_chromosomes = [chromosome.clone() for chromosome in genome.chromosomes]
    for chromosome in dup_chromosomes:
        genome.chromosomes.append(chromosome)


def apply_random_events(source_genome, dist_param, current_insertion_gene):
    """

    :param dist_param: DistanceParameters
    :param source_genome: Genome
    :param current_insertion_gene: Int
    :return:
    """
    insertion_length_range = range(1, dist_param.indel_length + 1)
    deletion_length_range = range(1, dist_param.indel_length + 1)
    duplication_length_range = range(1, dist_param.duplication_length + 1)

    # choose events and apply:
    event_count = {event: 0 for event in EventType.all}
    events = np.random.choice(
        [EventType.REARRANGEMENT, EventType.INSERTION, EventType.DELETION, EventType.DUPLICATION],
        dist_param.self.event_number,
        p=[dist_param.rearrangement_p, dist_param.insertion_p, dist_param.deletion_p, dist_param.duplication_p])

    for event in events:  # number of events, can be weighted by 'scaling' parameters
        if event == EventType.REARRANGEMENT:
            apply_random_rearrangement(param=dist_param, genome=source_genome)
        elif event == EventType.DELETION:
            apply_random_deletion(genome=source_genome, deletion_length_range=deletion_length_range)
        elif event == EventType.INSERTION:
            current_insertion_gene = apply_random_insertion(genome=source_genome,
                                                            gene=current_insertion_gene,
                                                            insertion_length_range=insertion_length_range)
        elif event == EventType.DUPLICATION:
            apply_random_segmental_duplication(genome=source_genome,
                                               duplication_length_range=duplication_length_range,
                                               interchromosome_p=dist_param.duplication_interchromosome_p)
        else:
            raise RuntimeError("Unknown evolutionary event.")
        event_count[event] += 1

    return event_count, current_insertion_gene
