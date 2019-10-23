import json
import os

LINEAR_END_CHR = "$"
CIRCULAR_END_CHR = "@"


def save_simulation(genomes, path_to_dir):
    if not os.path.exists(path_to_dir):
        os.makedirs(path_to_dir)

    write_genomes_to_file(genomes, os.path.join(path_to_dir, "blocks.txt"))
    for genome in genomes:
        write_genomes_to_file([genome], os.path.join(path_to_dir, genome.name + ".txt"))


def write_genomes_to_file(genomes, filename, write_chr_line=True):
    """
    Write genomes in a file with GRIMM format.
    """
    assert isinstance(genomes, list)

    with open(filename, "w") as f:
        for genome in genomes:
            f.write(">{0}\n".format(genome.name))
            for idx, chromosome in enumerate(genome.chromosomes):
                if write_chr_line:
                    f.write("# chr%d\n" % (idx + 1))
                f.write("%s %s\n" % (" ".join([str(gene) for gene in chromosome.gene_order]),
                                     CIRCULAR_END_CHR if chromosome.circular else LINEAR_END_CHR))


def write_simulation_parameters(param, output_path):
    with open(os.path.join(output_path, "params.cfg"), "w") as f:
        json.dump(param.__dict__, f, sort_keys=True, indent=4)
