# -*- coding:utf-8 -*-
import pandas as pd
import optparse
import yaml
import vcf
import sys
import os


def print_usage(option, opt, value, parser):
    usage_message = r"""
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
  Sensitivity and Specificity Statistics for Standards.
  python Standards_statistics.py -c StandardStat.yaml -s sampleID -g StandardsName \
  -v vcf_file -o outdir
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
    """
    print(usage_message)
    sys.exit()


def yaml_read(yaml_file):
    with open(yaml_file, 'r') as y:
        yaml_dic = yaml.load(y, Loader=yaml.FullLoader)
    return yaml_dic


def get_vcf_dp_gt_ratio_flt(vcf_file):
    chrom, pos, ref, alt, depth, genotype, ratio, flt, mutype = [], [], [], [], [], [], [], [], []
    vcf_reader = vcf.Reader(filename=vcf_file)
    for record in vcf_reader:
        try:
            if record.samples[0]['AD']:
                if record.samples[0]['GT'] != '1/2':
                    ratio.append(float(record.samples[0]['AD'][1])/(float(record.samples[0]['AD'][0])+float(record.samples[0]['AD'][1])))
                else:
                    ratio.append(float(record.samples[0]['AD'][2])/(float(record.samples[0]['AD'][1])+float(record.samples[0]['AD'][2])))
                chrom.append(record.CHROM)
                pos.append(record.POS)
                ref.append(record.REF)
                if len(record.ALT) == 1:
                    alt.append(str(record.ALT[0]))
                else:
                    alt.append(','.join([str(record.ALT[0]), str(record.ALT[1])]))
                depth.append(record.INFO['DP'])
                genotype.append(record.samples[0]['GT'])
                if record.FILTER:
                    flt.append('Filter')
                else:
                    flt.append('Pass')
                if record.is_snp:
                    mutype.append("snp")
                elif record.is_indel:
                    mutype.append("indel")
                else:
                    mutype.append("")
        except:
            pass
    dpgt_dict = {'#CHROM': chrom, 'POS': pos, 'REF': ref, 'ALT': alt, 'DP': depth, 'GT': genotype, 'Ratio': ratio, 'FILTER': flt, 'Type': mutype}
    dpgt_df = pd.DataFrame(dpgt_dict, columns=['#CHROM', 'POS', 'REF', 'ALT', 'DP', 'GT', 'Ratio', 'FILTER', 'Type'])
    return dpgt_df


def change_vcf_gt(vcf_dpgt):
    vcf_dpgt.loc[vcf_dpgt['GT'] == '0|1', 'GT'] = '0/1'
    vcf_dpgt.loc[vcf_dpgt['GT'] == '1|0', 'GT'] = '0/1'
    vcf_dpgt.loc[vcf_dpgt['GT'] == '1|1', 'GT'] = '1/1'
    vcf_dpgt.loc[vcf_dpgt['GT'] == '1|2', 'GT'] = '1/2'
    vcf_dpgt.loc[vcf_dpgt['GT'] == '2|1', 'GT'] = '1/2'
    return vcf_dpgt


def gold_sensitivity(snp_het, snp_hom, indel_het, indel_hom, goldid, sample_df, sampleid):
    sta = {}
    cols = ["#CHROM", "POS", "REF", "ALT", "GT"]
    gold_snp_het = pd.read_csv(snp_het, sep='\t')
    gold_snp_hom = pd.read_csv(snp_hom, sep='\t')
    gold_indel_het = pd.read_csv(indel_het, sep='\t')
    gold_indel_hom = pd.read_csv(indel_hom, sep='\t')
    sta['gold_snp_het'] = gold_snp_het.shape[0]
    sta['gold_snp_hom'] = gold_snp_hom.shape[0]
    sta['gold_indel_het'] = gold_indel_het.shape[0]
    sta['gold_indel_hom'] = gold_indel_hom.shape[0]
    sample_df = change_vcf_gt(sample_df)
    sample_df_snp_het = sample_df[(sample_df['Type'] == 'snp') & (sample_df['GT'] != '1/1')].copy()
    sample_df_snp_hom = sample_df[(sample_df['Type'] == 'snp') & (sample_df['GT'] == '1/1')].copy()
    sample_df_indel_het = sample_df[(sample_df['Type'] == 'indel') & (sample_df['GT'] != '1/1')].copy()
    sample_df_indel_hom = sample_df[(sample_df['Type'] == 'indel') & (sample_df['GT'] == '1/1')].copy()
    sta[goldid + '_snp_het'] = sample_df_snp_het.shape[0]
    sta[goldid + '_snp_hom'] = sample_df_snp_hom.shape[0]
    sta[goldid + '_indel_het'] = sample_df_indel_het.shape[0]
    sta[goldid + '_indel_hom'] = sample_df_indel_hom.shape[0]
    sample_df_snp_het_hit = pd.merge(gold_snp_het, sample_df_snp_het, on=cols, how='inner')
    sample_df_snp_hom_hit = pd.merge(gold_snp_hom, sample_df_snp_hom, on=cols, how='inner')
    sample_df_indel_het_hit = pd.merge(gold_indel_het, sample_df_indel_het, on=cols, how='inner')
    sample_df_indel_hom_hit = pd.merge(gold_indel_hom, sample_df_indel_hom, on=cols, how='inner')
    sta[goldid + '_snp_het_hit'] = sample_df_snp_het_hit.shape[0]
    sta[goldid + '_snp_hom_hit'] = sample_df_snp_hom_hit.shape[0]
    sta[goldid + '_indel_het_hit'] = sample_df_indel_het_hit.shape[0]
    sta[goldid + '_indel_hom_hit'] = sample_df_indel_hom_hit.shape[0]
    sta['sample'] = sampleid
    sta_cols = ['sample', 'gold_snp_het', 'gold_snp_hom', 'gold_indel_het', 'gold_indel_hom', goldid + '_snp_het', goldid + '_snp_hom',
                goldid + '_indel_het', goldid + '_indel_hom', goldid + '_snp_het_hit', goldid + '_snp_hom_hit',
                goldid + '_indel_het_hit', goldid + '_indel_hom_hit']
    sta_df = pd.DataFrame(sta, columns=sta_cols, index=[0])
    return sta_df


def gold_specificity(tn_bed, fp_vcf):
    tn_sta = {}
    tn_bed_df = pd.read_csv(tn_bed, sep='\t', header=None)
    tn_bed_df.columns = ['#Chr', 'Start', 'Stop']
    tn_bed_df['size'] = tn_bed_df['Stop'] - tn_bed_df['Start'] + 1
    tn_sta['TN_bed_size'] = tn_bed_df['size'].sum()
    try:
        tmp_vcf = pd.read_csv(fp_vcf, header=None, sep='\t')
        tn_sta['var_in_TN_bed'] = tmp_vcf.shape[0]
        tn_sta['snp_var_in_TN_bed'] = tmp_vcf[(tmp_vcf[3].map(len) == 1) & (tmp_vcf[4].map(len) == 1)].shape[0]
        tn_sta['indel_var_in_TN_bed'] = tmp_vcf[(tmp_vcf[3].map(len) != 1) | (tmp_vcf[4].map(len) != 1)].shape[0]
    except:
        tn_sta['var_in_TN_bed'] = 0
        tn_sta['snp_var_in_TN_bed'] = 0
        tn_sta['indel_var_in_TN_bed'] = 0
    tn_sta['TN'] = tn_sta['TN_bed_size'] - tn_sta['var_in_TN_bed']
    sta_cols = ['TN_bed_size', 'TN', 'var_in_TN_bed', 'snp_var_in_TN_bed', 'indel_var_in_TN_bed']
    sta_df = pd.DataFrame(tn_sta, columns=sta_cols, index=[0])
    return sta_df


def gold_stat(sen_df, spe_df, goldid):
    sta = sen_df.join(spe_df)
    sta['snp_het_sensitivity'] = sta[goldid + '_snp_het_hit'] / sta['gold_snp_het']
    sta['snp_hom_sensitivity'] = sta[goldid + '_snp_hom_hit'] / sta['gold_snp_hom']
    sta['indel_het_sensitivity'] = sta[goldid + '_indel_het_hit'] / sta['gold_indel_het']
    sta['indel_hom_sensitivity'] = sta[goldid + '_indel_hom_hit'] / sta['gold_indel_hom']
    sta['snp_sensitivity'] = (sta[goldid + '_snp_het_hit'] + sta[goldid + '_snp_hom_hit']) / (
            sta['gold_snp_het'] + sta['gold_snp_hom'])
    sta['indel_sensitivity'] = (sta[goldid + '_indel_het_hit'] + sta[goldid + '_indel_hom_hit']) / (
            sta['gold_indel_het'] + sta['gold_indel_hom'])
    sta['snp_TN'] = sta['TN_bed_size'] - sta['snp_var_in_TN_bed']
    sta['indel_TN'] = sta['TN_bed_size'] - sta['indel_var_in_TN_bed']
    sta['snp_specificity'] = sta['snp_TN'] / sta['TN_bed_size']
    sta['indel_specificity'] = sta['indel_TN'] / sta['TN_bed_size']
    sta['specificity'] = sta['TN'] / sta['TN_bed_size']
    sta['FN'] = sta['gold_snp_het'] - sta[goldid + '_snp_het_hit'] + sta['gold_snp_hom'] - sta[goldid + '_snp_hom_hit'] + \
                sta['gold_indel_het'] - sta[goldid + '_indel_het_hit'] + sta['gold_indel_hom'] - sta[goldid + '_indel_hom_hit']
    sta['TP'] = sta[goldid + '_snp_het_hit'] + sta[goldid + '_snp_hom_hit'] + sta[goldid + '_indel_het_hit'] + \
                sta[goldid + '_indel_hom_hit']
    sta['FP'] = sta['var_in_TN_bed']
    sta['PPV'] = sta['TP'] / (sta['TP'] + sta['FP'])
    sta['Accuracy'] = (sta['TN'] + sta['TP']) / (sta['TN'] + sta['TP'] + sta['FN'] + sta['FP'])
    sta['snp_FN'] = sta['gold_snp_het'] - sta[goldid + '_snp_het_hit'] + sta['gold_snp_hom'] - sta[goldid + '_snp_hom_hit']
    sta['indel_FN'] = sta['gold_indel_het'] - sta[goldid + '_indel_het_hit'] + sta['gold_indel_hom'] - sta[goldid + '_indel_hom_hit']
    sta['snp_TP'] = sta[goldid + '_snp_het_hit'] + sta[goldid + '_snp_hom_hit']
    sta['indel_TP'] = sta[goldid + '_indel_het_hit'] + sta[goldid + '_indel_hom_hit']
    sta['snp_FP'] = sta['snp_var_in_TN_bed']
    sta['indel_FP'] = sta['indel_var_in_TN_bed']
    sta['snp_accuracy'] = (sta['snp_TN'] + sta['snp_TP']) / (sta['snp_TN'] + sta['snp_TP'] + sta['snp_FN'] + sta['snp_FP'])
    sta['indel_accuracy'] = (sta['indel_TN'] + sta['indel_TP']) / (sta['indel_TN'] + sta['indel_TP'] + sta['indel_FN'] + sta['indel_FP'])
    return sta


if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('-u', '--usage', help='print more info on how to use this script', action="callback", callback=print_usage)
    parser.add_option('-c', '--config', dest='config', help='config file', default=None, metavar='file')
    parser.add_option('-s', '--sample', dest='sample', help='sample ID', default=None, metavar='string')
    parser.add_option('-g', '--gold', dest='gold', help='Name of standards', default=None, metavar='string')
    parser.add_option('-v', '--vcf', dest='vcf_file', help='vcf file', default=None, metavar='file')
    parser.add_option('-o', '--outdir', dest='outdir', help='work dir', default=None, metavar='string')
    (opts, args) = parser.parse_args()
    config = opts.config
    sample = opts.sample
    gold = opts.gold
    vcf_file = opts.vcf_file
    outdir = opts.outdir
    config_dic = yaml_read(config)
    df_vcf = get_vcf_dp_gt_ratio_flt(vcf_file)
    df_vcf.to_csv(outdir + "/" + sample + ".dpgt", sep='\t', index=False)
    com = config_dic['tabix'] + " -R " + config_dic['tn_bed'] + " " + vcf_file + ">" + os.path.join(outdir, sample + ".fp.vcf")
    os.system(com)
    df_sen = gold_sensitivity(config_dic['snp_het'], config_dic['snp_hom'], config_dic['indel_het'], config_dic['indel_hom'], gold, df_vcf, sample)
    df_spe = gold_specificity(config_dic['tn_bed'], os.path.join(outdir, sample + ".fp.vcf"))
    stat = gold_stat(df_sen, df_spe, gold)
    stat.T.to_excel(outdir + "/" + sample + ".sensitivity.xlsx", header=None)
