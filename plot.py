import pandas as pd
import matplotlib.pyplot as plt
import common as cmn
import os
import logging


def show_plot(file_path, gene_a=cmn.PTEN_GENE_CODE, gene_b=cmn.PTEN_GENE_CODE):
    df = pd.read_csv(file_path, sep="\t", index_col=0)
    df.T.plot(y=[gene_a, gene_b], style=".")
    plt.show()


def corr_gene_pair(dt, gene_a=cmn.PTEN_GENE_CODE, gene_b=cmn.PTEN_GENE_CODE):
    sd = dt[[gene_a, gene_b]]
    return sd.corr().values[0][1]


def corr_genes(r, g):
    """

    :return:
    """
    norm = "normalized"
    if r:
        norm = "raw"

    for cancer_type in cmn.list_dir(cmn.DL_DIR):
        for gender in cmn.list_dir(os.path.join(cmn.DL_DIR, cancer_type)):
            logging.info("Creating correlation files for '%s' > '%s' ..." % (cancer_type, gender))

            samples_path = os.path.join(cmn.DL_DIR, cancer_type, gender)

            integrated_file = cmn.INTEGRATED_FNAME % (cancer_type.lower(), gender.lower(), "rna", norm)
            integrated_file_path = os.path.join(samples_path, integrated_file)

            corr_file = cmn.CORR_FNAME % (cancer_type.lower(), gender.lower(), "rna", norm)
            corr_file_path = os.path.join(samples_path, corr_file)

            if os.path.isfile(integrated_file_path):
                logging.info("Found integrated file, processing '%s' ...\n" % integrated_file)

                df = pd.read_csv(integrated_file_path, sep="\t", index_col=0)
                dt = df.T
                gene_names = df.index
                num_genes = len(gene_names)
                n = 0

                if g == "all" or g == "" or not g.startswith("E"):
                    # TODO
                    logging.info("TODO: Creating correlation index for all genes in this file ...")
                else:
                    logging.info("Creating correlation index for gene [%s] with all other genes in this file ..." % g)

                    out_lines = ["0\t%s" % g]
                    for gene in gene_names:
                        logging.info("Processing gene [%s (%s out of %s), in '%s' > '%s'] ..."
                                     % (gene, n, num_genes, cancer_type, gender))
                        out_lines.append("%s\t%s" % (gene, corr_gene_pair(dt, g, gene)))
                        n += 1

                    logging.info("Writing values to file '%s' ..." % corr_file)
                    with open(corr_file_path, 'w') as out_file:
                        out_file.write("\n".join(out_lines))

                    logging.info("Done.\n")

            else:
                logging.error("No integrated (%s) RNA file found for '%s' > '%s'.\n"
                              % (norm, cancer_type, gender))



def run(r=False, g=cmn.PTEN_GENE_CODE):
    """
    run this module

    :return:
    """
    corr_genes(r, g)
