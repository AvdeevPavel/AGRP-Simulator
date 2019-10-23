#!/usr/bin/env python3

import argparse
import os
import shutil


def generate_double_distance_problems(config, out_dir):
    pass


def generate_guided_halving_problems(config, out_dir):
    pass


def generate_median_problems(config, out_dir):
    pass


def get_immediate_subdirectories(a_dir):
    """
    This function get subdirectories
    """
    return ((os.path.join(a_dir, name), name) for name in os.listdir(a_dir) if os.path.isdir(os.path.join(a_dir, name)))


def get_immediate_files(a_dir):
    """
    This function get files
    """
    return ((os.path.join(a_dir, name), name) for name in os.listdir(a_dir) if
            os.path.isfile(os.path.join(a_dir, name)))


def dojob_split(name_directory):
    # logging.info('Let us do the ILP')

    for path, name in get_immediate_subdirectories(name_directory):
        for subpath, subname in get_immediate_subdirectories(path):
            if os.path.isfile(subpath):
                continue

            for name_dir in ["heur_ilp"]: # "force_ilp", "rgghp_ilp", "heur_ilp"]:
                force_ilp_dir = os.path.join(subpath, name_dir)
                if not os.path.exists(force_ilp_dir):
                    os.makedirs(force_ilp_dir)

                # rgghp_ilp_dir = os.path.join(subpath, "rgghp_ilp")
                # if not os.path.exists(rgghp_ilp_dir):
                #     os.makedirs(rgghp_ilp_dir)
                #
                # heur_ilp_dir = os.path.join(subpath, "heur_ilp")
                # if not os.path.exists(heur_ilp_dir):
                #     os.makedirs(heur_ilp_dir)

                shutil.copyfile(os.path.join(subpath, "A.txt"), os.path.join(force_ilp_dir, "A.gen"))
                shutil.copyfile(os.path.join(subpath, "B.txt"), os.path.join(force_ilp_dir, "B.gen"))
                shutil.copyfile(os.path.join(subpath, "R.txt"), os.path.join(force_ilp_dir, "R.gen"))
                shutil.copyfile(os.path.join(subpath, "params.cfg"), os.path.join(force_ilp_dir, "params.cfg"))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Simulates rearrangement evolution w/ a WGD event")
    parser.add_argument("-wgd", "--WGD_type", type=str, default=None,
                        help="Type of WGD event, either 'DD' for double distance or 'GH' for genome halving or None. ")
    parser.add_argument("-n", "--num_genes", type=int, default=100, help="Number of genes in the root genome.")
    parser.add_argument("-c", "--num_chr", type=int, default=5, help="Number of chromosomes in the root genome.")
    parser.add_argument("-ev", "--num_ev", type=int, default=10, help="Number of rearrangement events before WGD.")
    parser.add_argument("-o", "--output", type=str, default="sim", help="Name of the output folder.")
    parser.add_argument("-dp", "--deletion_p", type=float, default=0.0,
                        help="Percentage of deletions before WGD, from 0 to 1.0")
    parser.add_argument("-ip", "--insertion_p", type=float, default=0.0,
                        help="Percentage of insertions before WGD, from 0 to 1.0")
    parser.add_argument("-dl", "--duplication_length", type=int, default=5,
                        help="Maximum length of duplication events pre WGD.")
    parser.add_argument("-dup_p", "--duplication_p", type=float, default=0.0,
                        help="Percentage of duplications, from 0 to 1.0 pre WGD")
    parser.add_argument("-dl2", "--duplication_length2", type=int, default=5,
                        help="Maximum length of duplication event post WGD.")
    parser.add_argument("-dup_p2", "--duplication_p2", type=float, default=0.0,
                        help="Percentage of duplications, from 0 to 1.0 post WGD")
    parser.add_argument("-il", "--indel_length", type=int, default=5,
                        help="Maximum size of indel event in genes (pre-WGD).")
    parser.add_argument("-ev2", "--num_ev2", type=int, default=10, help="Number of rearrangement events after WGD.")
    parser.add_argument("-dp2", "--deletion_p2", type=float, default=0.0,
                        help="Percentage of deletions after WGD, from 0 to 1.0")
    parser.add_argument("-ip2", "--insertion_p2", type=float, default=0.0,
                        help="Percentage of insertions after WGD, from 0 to 1.0")
    parser.add_argument("-il2", "--indel_length2", type=int, default=5,
                        help="Maximum size of indel event in genes after WGD.")
    parser.add_argument("-idup_p", "--inter_duplication_p", type=float, default=0.0,
                        help="Percentage of segmental duplications that occur across two chromosomes")

    param = parser.parse_args()

