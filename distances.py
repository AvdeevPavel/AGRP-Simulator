from events import apply_wgd, apply_random_events


class SimGenomeParameters:
    def __init__(self, num_genes=0, num_chr=0):
        self.num_genes = num_genes
        self.num_chr = num_chr


def simulate_regular_distance(source_genome, name, dist_parameters, current_insertion_gene):
    target_genome = source_genome.clone(name=name)
    event_count, current_insertion_gene = apply_random_events(dist_param=dist_parameters, source_genome=target_genome,
                                                              current_insertion_gene=current_insertion_gene)
    return target_genome, event_count, current_insertion_gene


def simulate_distance_with_wgd(source_genome, name, dist_parameters, current_insertion_gene):
    target_genome = source_genome.clone(name=name)
    apply_wgd(genome=target_genome)
    event_count, current_insertion_gene = apply_random_events(dist_param=dist_parameters, source_genome=target_genome,
                                                              current_insertion_gene=current_insertion_gene)
    return target_genome, event_count, current_insertion_gene
