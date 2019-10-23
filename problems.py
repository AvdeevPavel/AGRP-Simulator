from distances import simulate_distance_with_wgd, simulate_regular_distance
from genome import Genome
from writer import save_simulation


def simulate_median_instance(genome_parameters, dist_parameters, path_to_dir):
    current_insertion_gene = genome_parameters.num_genes + 1
    median_genome = Genome.identity(genome_parameters.num_genes, genome_parameters.num_chr, name="median")

    genomes = []
    for i in range(3):
        genome, _, current_ins_gene = simulate_regular_distance(source_genome=median_genome,
                                                                dist_parameters=dist_parameters,
                                                                current_insertion_gene=current_insertion_gene,
                                                                name="P" + str(i + 1))
        genomes.append(genome)
        current_insertion_gene = current_ins_gene

    genomes.append(median_genome)

    save_simulation(genomes, path_to_dir)


def simulate_double_distance_instance(genome_parameters, dist_parameters, path_to_dir):
    current_insertion_gene = genome_parameters.num_genes + 1
    predup_genome = Genome.identity(genome_parameters.num_genes, genome_parameters.num_chr, name="R")
    a_genome, _, _ = simulate_distance_with_wgd(source_genome=predup_genome, dist_parameters=dist_parameters,
                                                current_insertion_gene=current_insertion_gene,
                                                name="A")

    save_simulation([predup_genome, a_genome], path_to_dir)


def simulate_guided_halving_instance(genome_parameters, dist_parameters_before, dist_parameters_after, path_to_dir):
    current_insertion_gene = genome_parameters.num_genes + 1
    b_genome = Genome.identity(genome_parameters.num_genes, genome_parameters.num_chr, name="B")

    predup_genome, _, current_ins_gene = simulate_regular_distance(source_genome=b_genome,
                                                                   dist_parameters=dist_parameters_before,
                                                                   current_insertion_gene=current_insertion_gene,
                                                                   name="R")

    a_genome, _, _ = simulate_distance_with_wgd(source_genome=predup_genome,
                                                dist_parameters=dist_parameters_after,
                                                current_insertion_gene=current_insertion_gene,
                                                name="A")

    save_simulation([b_genome, predup_genome, a_genome], path_to_dir)
