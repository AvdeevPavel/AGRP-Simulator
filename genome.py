from collections import Counter


class Genome:
    def __init__(self, name="", chr_list=None):
        if not chr_list:
            chr_list = []

        self.name = name
        self.chromosomes = list()

        for gene_order in chr_list:
            if isinstance(gene_order, Chromosome):
                self.chromosomes.append(gene_order)
            elif isinstance(gene_order, tuple):
                self.chromosomes.append(Chromosome(list(gene_order), circular=True))
            else:
                self.chromosomes.append(Chromosome(gene_order, circular=False))

    @staticmethod
    def identity(num_genes, num_chromosomes, name="Identity", circular=False):
        genes_per_chr = num_genes / num_chromosomes
        return Genome(name, [Chromosome(range(c * genes_per_chr + 1, (c + 1) * genes_per_chr + 1), circular=circular)
                             for c in range(num_chromosomes)])

    def add_chromosome(self, chromosome):
        self.chromosomes.append(chromosome)

    def n_chromosomes(self):
        return len(self.chromosomes)

    def gene_count(self):
        return Counter([abs(x) for chromosome in self.chromosomes for x in chromosome.gene_order])

    def gene_set(self):
        return set([abs(x) for chromosome in self.chromosomes for x in chromosome.gene_order])

    def adjacency_set(self):
        return set.union(*[chromosome.adjacency_set() for chromosome in self.chromosomes])

    def clone(self, name=None):
        g = Genome(name=name if name is not None else self.name)
        for chromosome in self.chromosomes:
            g.add_chromosome(chromosome.clone())
        return g

    def common_adjacencies(self, g):
        return self.adjacency_set().intersection(g.adjacency_set())

    def __str__(self):
        return ">" + self.name + "\n" + "[" + ", ".join([str(c) for c in self.chromosomes]) + "]"

    def build_extremity_order_lists(self):
        ext_order_list = []
        for idx, chrom in enumerate(self.chromosomes):
            ext_order_list.append(chrom.build_chromosome_ext_order())
        return ext_order_list


class Chromosome:
    def __init__(self, gene_order, circular=False):
        self.circular = circular
        self.gene_order = gene_order

    def __iter__(self):
        return self.gene_order.__iter__()

    def next(self):
        return self.gene_order.next()

    def __str__(self):
        delimiter = ("(", ")") if self.circular else ("[", "]")
        return "%s%s%s" % (delimiter[0], ", ".join([str(x) for x in self.gene_order]), delimiter[1])

    def length(self):
        return len(self.gene_order)

    def adjacency_set(self):
        def adjacency(g1, g2):
            ext1 = 2 * abs(g1) - 1 if g1 < 0 else 2 * abs(g1)
            ext2 = 2 * abs(g2) if g2 < 0 else 2 * abs(g2) - 1
            return (ext1, ext2) if ext1 < ext2 else (ext2, ext1)

        s = {adjacency(gene1, gene2) for gene1, gene2 in zip(self.gene_order, self.gene_order[1:])}
        if self.circular:
            s.add(adjacency(self.gene_order[-1], self.gene_order[0]))
        return s

    def clone(self):
        return Chromosome(list(self.gene_order), self.circular)
