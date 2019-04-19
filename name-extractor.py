import pandas as pd
import numpy as np
import argparse as arg
import time


def match_rows(haystack, h_search_col, needle_search_col):
	return haystack[haystack[h_search_col].isin(needle_search_col)]


def combine_attr(df, att_col, delim=";"):
	attributes = split_attribute(df[att_col], delim)
	new_df = df.drop(df.columns[-1], axis=1)
	combined = pd.concat([new_df, attributes], axis=1)
	return combined


def split_attribute(attribute, delim):
	split_cols = attribute.str.split(delim, expand=True)
	split_cols = split_cols.drop(split_cols.columns[-1], axis=1)
	col_name_row = split_cols.iloc[0].str.strip()
	col_names = col_name_row.str.extract(r"(\w*)")[0].tolist()
	split_cols.columns = col_names

	for col in split_cols:
		split_cols[col] = split_cols[col].str.extract(r'\"(.*)\"')

	return split_cols


def unmatched_rows(matched_search_col, needle_search_col):
	combined = pd.concat([matched_search_col, needle_search_col])
	duplicated = combined.duplicated(keep=False)
	unmatched_indices = duplicated[duplicated == False].index.values
	return needle_search_col[unmatched_indices]


def get_filetype(fname):
	start = fname.rfind(".") + 1
	return fname[start:]


def is_num(n):
	try:
		int(n)
	except ValueError:
		return False
	return True


def extract_main(n_path, h_path, n_delim, h_delim, n_header, h_header, n_col, h_col, h_split_col, h_split_delim, out_path):
	n_ftype = get_filetype(n_path)
	h_ftype = get_filetype(h_path)

	print("Reading input files...")

	if n_ftype == "xlsx":
		needle = pd.read_excel(n_path)
	else:
		needle = pd.read_csv(n_path, delimiter=n_delim, header=n_header)

	if h_ftype == "xlsx":
		haystack = pd.read_excel(h_path)
	else:
		haystack = pd.read_csv(h_path, delimiter=h_delim, header=h_header)

	print("Done!")

	if h_split_col != None:
		print("Splitting column...")
		haystack = combine_attr(haystack, h_split_col, h_split_delim)
		print("Done!")

	print("Extracting rows...")
	extracted = match_rows(haystack, h_col, needle[n_col])
	print("Done!")

	extracted.to_csv(out_path, index=False, sep=h_delim, header=h_header)
	print("Extracted rows written to: \"{}\"".format(out_path))

	unmatched = unmatched_rows(extracted[h_col], needle[n_col])
	print("The following entries could not be extraced ({}/{}):".format(len(unmatched), len(needle)))
	for entry in unmatched.values:
		print(entry)


def cli_args():
	parser = arg.ArgumentParser(description="Extract rows of a delimited haystack file based on a search column from a needle file.")
	parser.add_argument("-ih", "--input_haystack", dest="h_in", required=True, help="Name of the input haystack file")
	parser.add_argument("-in", "--input_needle", dest="n_in", required=True, help="Name of the input needle file")
	parser.add_argument("-hc", "--haystack_col", dest="h_col", required=True, help="Column name to be searched in the haystack file")
	parser.add_argument("-nc", "--needle_col", dest="n_col", required=True, help="Column name of the needle file used to query the haystack file")
	parser.add_argument("-hd", "--haystack_delim", dest="h_delim", default=",", help="[OPTIONAL] Delimiter for the haystack file (default=\",\")")
	parser.add_argument("-nd", "--needle_delim", dest="n_delim", default=",", help="[OPTIONAL] Delimiter for the needle file (default=\",\")")
	parser.add_argument("-hh", "--haystack_header", dest="h_header", default=None, action="store_const", const="0", help="[OPTIONAL] Flag to indicate the presence of column headers in the haystack file")
	parser.add_argument("-nh", "--needle_header", dest="n_header", default=None, action="store_const", const="0", help="[OPTIONAL] Flag to indicate the presence of column headers in the needle file (default=\"None\")")
	parser.add_argument("-sc", "--split_col", dest="split_col", default=None, help="[OPTIONAL] Name of column to split in the haystack file (useful for GTF attribute column, default=\"None\")")
	parser.add_argument("-sd", "--split_delim", dest="split_delim", default=None, help="[OPTIONAL] Delimiter used to split column if --split_col is specified (default=\"None\")")
	parser.add_argument("-o", "--out_path", dest="out_path", default="extracted.csv", help="[OPTIONAL] Name of the output file (default=\"extracted.csv\")")
	args = parser.parse_args()

	if args.n_delim == "tab":
		args.n_delim = "\t"
	if args.h_delim == "tab":
		args.h_delim = "\t"

	if is_num(args.n_col):
		args.n_col = int(args.n_col)
	if is_num(args.h_col):
		args.h_col = int(args.h_col)
	if is_num(args.split_col):
		args.split_col = int(args.split_col)

	main_args = {
		"n_path": args.n_in,
		"h_path": args.h_in,
		"n_delim": args.n_delim,
		"h_delim": args.h_delim,
		"n_header": args.n_header,
		"h_header": args.h_header,
		"n_col": args.n_col,
		"h_col": args.h_col,
		"h_split_col": args.split_col,
		"h_split_delim": args.split_delim,
		"out_path": args.out_path
	}

	extract_main(**main_args)


cli_args()


"""
delim = "\t"
needle = pd.read_excel("C:\\Users\\dgall\\OneDrive\\Documents\\Research\\loci_names_master.xlsx")
haystack = pd.read_csv("correct_cuffmerge.gtf", header=None, delimiter="\t", nrows=1)
#needle = pd.read_csv("cuff_needle_test.csv")
#haystack = pd.read_csv("cuff_haystack_test.csv", header=None)
print("Haystack length:", len(haystack))

split_start = time.time()
print("Splitting GTF attributes...")
haystack = combine_attr(haystack, int("8"))
split_delta = time.time() - split_start
print("Done! Split execution time: %.5fs" % split_delta)

print()

m_start = time.time()
print("Matching needles to haystack file...")
matched = match_rows(haystack, "gene_id", needle['XLOC_id'])
m_delta = time.time() - m_start
print("Done! Matching execution time: %.5fs" % m_delta)

print()

#matched.to_csv("matched_out.csv", index=False)
print("Matches written to \"matched_out.csv\"")

print()

print("Matches to the following needles could not be found:")
unmatched = unmatched_rows(matched['gene_id'], needle['XLOC_id'])
print("Number of unmatched needles: {} (out of {})".format(len(unmatched), len(needle['XLOC_id'])))
for entry in unmatched.values:
	print(entry)
"""