# -*- coding:utf-8 -*-
import pandas as pd
import optparse
import sys


def print_usage(option, opt, value, parser):
    usage_message = r"""
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
  Generate Standards Set
  python Standards_set.py -p /path/to/work -v Standards.vcf.gz -g StandardsName \
  --col_name GT_info_column --skip_row Standards_vcf_header_rows
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
    """
    print(usage_message)
    sys.exit()


def get_vcf_info(vcf, skips, col):
    vcf_df = pd.read_csv(vcf, skiprows=range(skips), sep='\t', dtype={'#CHROM': str}, compression='gzip')
    vcf_df['#CHROM'] = 'chr' + vcf_df['#CHROM']
    vcf_df['GT'] = vcf_df[col].str.split(':').str[0]
    vcf_df['GT'] = vcf_df['GT'].str.replace('|', '/')
    vcf_df.loc[vcf_df['GT'] == '1/0', 'GT'] = '0/1'
    vcf_df.loc[vcf_df['GT'] == '2/1', 'GT'] = '1/2'
    vcf_df['Type'] = 'indel'
    vcf_df.loc[(vcf_df['GT'] != '1/2') & (vcf_df['REF'].str.len() == 1) & (vcf_df['ALT'].str.len() == 1), 'Type'] = 'snv'
    vcf_df.loc[(vcf_df['GT'] == '1/2') & (vcf_df['REF'].str.len() == 1) & (vcf_df['ALT'].str.len() == 3), 'Type'] = 'snv'
    return vcf_df[["#CHROM", "POS", "ID", "REF", "ALT", "GT", "Type"]]


def split_hom_het(vcf_df):
    vcf_hom = vcf_df[vcf_df['GT'] == '1/1'].copy()
    vcf_het = vcf_df[vcf_df['GT'] != '1/1'].copy()
    return vcf_hom, vcf_het


def write_vcf_info(vcf_df, wdir, sample, zy):
    vcf_df_snv = vcf_df[vcf_df['Type'] == 'snv'].copy()
    del vcf_df_snv['Type']
    vcf_df_indel = vcf_df[vcf_df['Type'] == 'indel'].copy()
    del vcf_df_indel['Type']
    vcf_df_snv.to_csv(wdir + "/" + sample + ".snv." + zy + '.txt', sep='\t', index=False)
    vcf_df_indel.to_csv(wdir + "/" + sample + ".indel." + zy + '.txt', sep='\t', index=False)


if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('-u', '--usage', help='print more info on how to use this script', action="callback", callback=print_usage)
    parser.add_option('-v', '--vcf', dest='vcf_file', help="Standards vcf file", default=None, type='string')
    parser.add_option('--skip_row', dest='skip_row', help="Standards vcf file header rows", type=int)
    parser.add_option('--col_name', dest='col_name', help="Standards vcf file GT information column", type='string')
    parser.add_option('-g', '--gold', dest='gold', help='Name of standards', default=None, type='string')
    parser.add_option('-p', '--pwd', dest='pwd', help="work dir", default=None, type='string')
    (opts, args) = parser.parse_args()
    vcf_file = opts.vcf_file
    skip_row = opts.skip_row
    out_file = opts.out_file
    col_name = opts.col_name
    pwd = opts.pwd
    vcf_bed = get_vcf_info(vcf_file, skip_row, col_name)
    vcf_bed_hom, vcf_bed_het = split_hom_het(vcf_bed)
    write_vcf_info(vcf_bed_hom, pwd, out_file, 'hom')
    write_vcf_info(vcf_bed_het, pwd, out_file, 'het')
