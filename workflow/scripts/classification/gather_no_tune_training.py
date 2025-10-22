#!/usr/bin/env python

import csv

def combine_tsv_files(file_list, outfile, add_col, missing_col):
    with open(outfile, 'w', newline='') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')

        for input_file in file_list:
            with open(input_file, 'r') as in_file:
                tsv_reader = csv.reader(in_file, delimiter='\t')
                
                header = next(tsv_reader)
                header_missing_col = -1

                if  add_col is not None:
                    header.append("primary_model")

                if missing_col is not None:
                    if missing_col not in header:
                        missing_col_index = len(header)
                    else:
                        missing_col_index = header.index(missing_col)
                        header_missing_col = header.pop(missing_col_index)
                    header.append(missing_col)

                if file_list.index(input_file) == 0:
                    tsv_writer.writerow(header)

                if add_col is not None:
                    if missing_col is not None:
                        if header_missing_col == -1:
                            for row in tsv_reader:
                                row.append(add_col)
                                row.append("FALSE")
                                tsv_writer.writerow(row)
                        else:
                            for row in tsv_reader:
                                missing_col_content = row.pop(missing_col_index)
                                row.append(add_col)
                                row.append(missing_col_content)
                                tsv_writer.writerow(row)
                    else:
                        for row in tsv_reader:
                            row.append(add_col)
                            tsv_writer.writerow(row)
                else:
                    if missing_col is not None:
                        if header_missing_col == -1:
                            for row in tsv_reader:
                                row.append("FALSE")
                                tsv_writer.writerow(row)
                        else:
                            for row in tsv_reader:
                                missing_col_content = row.pop(missing_col_index)
                                row.append(missing_col_content)
                                tsv_writer.writerow(row)
                    else:
                        for row in tsv_reader:
                            tsv_writer.writerow(row)

def main():
    prim_mod = snakemake.params.get("prim_mod", None)
    ## There is a missing column in some files. If this column is specified, fill in FALSE if missing
    fix_missing_col = snakemake.params.get("missing_col", None)

    combine_tsv_files(snakemake.input["probs"], snakemake.output["gathered_probs"], prim_mod, fix_missing_col)
    combine_tsv_files(snakemake.input["assessment_probs"], snakemake.output["gathered_assessment"], prim_mod, fix_missing_col)
    combine_tsv_files(snakemake.input["model_specs"], snakemake.output["gathered_model_specs"], prim_mod, fix_missing_col)
    combine_tsv_files(snakemake.input["feature_importance"], snakemake.output["gathered_fi"], prim_mod, fix_missing_col)


if __name__ == "__main__":
    main()
