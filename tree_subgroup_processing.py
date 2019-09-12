from ete3 import Tree
from Bio import SeqIO
import vcfpy
from argparse import (ArgumentParser, FileType)


def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(
        description='BIO_Hansel Scheme developement')

    parser.add_argument('--in_vcf', type=str, required=True, help='VCF File of SNPs')
    parser.add_argument('--in_nwk', type=str, required=True, help='Newick Tree of strains')
    parser.add_argument('--reference', type=str, required=True, help='Reference fasta sequence from VCF')
    parser.add_argument('--outdir', type=str, required=False, help='Output Directory to put results')
    parser.add_argument('--min_snps',type=int, required=False, default=2,help="Number of cannonical SNPs required to define a group")
    parser.add_argument('--min_members', type=int, required=False, default=5,
                        help="Minimum number of members to define a group")
    return parser.parse_args()

def check_snp(in_group,memberships,position,tile_size,ref_seq,valid_groups,valid_ranks,mid_point,ref_base,alt_base,count_ref,count_alt):
    candidate_partitions = dict()

    sample1 = in_group[0]
    max_rank = len(memberships[sample1])
    #print(memberships)

    for i in range(0, max_rank):

        heir_membership = dict()
        for sample2 in in_group:

            #print(memberships[sample2][i])
            heir_membership[memberships[sample2][i]] = ''

            start = position - int(tile_size / 2) - 1
            end = position + int(tile_size / 2)
            positive_tile = ''
            negative_tile = ''
            #positive_tile = ref_seq[start:end]
            #negative_tile = positive_tile
            #positive_tile = "{}{}{}".format(positive_tile[0:int(tile_size / 2) - 1],alt_base,positive_tile[int(tile_size / 2):tile_size])
            #positive_tile = ''.join([positive_tile[0:int(tile_size / 2) - 1],alt_base,positive_tile[int(tile_size / 2):tile_size]])
            #positive_tile = list(positive_tile)
            #positive_tile[mid_point] = alt_base
            #positive_tile = ''.join(str(v) for v in positive_tile)


        #print(heir_membership)

        if 0 in heir_membership:
            heir_membership = dict()

        if len(heir_membership) == 1 and 0 not in heir_membership:
            group_id = memberships[sample1][i]
            all_group_members = list()
            for n in memberships:
                if memberships[n][i] == group_id:
                    all_group_members.append(n)

            #print("{}\t{}\t{}".format(i,len(in_group),len(all_group_members)))
            if len(set().union(all_group_members + in_group)) == len(in_group):
                valid_ranks.append(i)
                valid_groups.append(group_id)
                temp = negative_tile
                negative_tile = positive_tile
                positive_tile = temp




            if len(valid_ranks) > 0:
                for i in range(0, len(valid_ranks)):
                    if valid_ranks[i] not in candidate_partitions:
                        candidate_partitions[valid_ranks[i]] = list()


                    candidate_partitions[valid_ranks[i]].append({"position": position,
                                                                 "ref_base": ref_base,
                                                                 "alt_base": alt_base,
                                                                 "ref_count": count_ref,
                                                                 "alt_count": count_alt,
                                                                 "positive_tile": positive_tile,
                                                                 "negative_tile": negative_tile,
                                                                 "group_id": valid_groups[i],
                                                                 "rank_id":valid_ranks[i]})

                    #print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(position,ref_base,alt_base,count_ref,count_alt,','.join(str(v) for v in valid_ranks),','.join(str(v) for v in valid_groups),"".join(positive_tile),"".join(negative_tile)))
    #print(candidate_partitions)
    return candidate_partitions


def parse_tree(tree_file):
    # Load a tree structure from a newick file.
    t = Tree(tree_file)

    #Need this otherwise groups derived from the tree are inaccurate
    t.resolve_polytomy()
    #Need to force root for consistency but should modify this behavior to support defined root

    root = t.get_midpoint_outgroup()
    #t.set_outgroup(root)
    return t


def get_tree_groups(ete3_tree_obj):
    level_rank = 0
    memberships = dict()
    group = 1

    for node in ete3_tree_obj.iter_descendants("levelorder"):
        names = node.get_leaf_names()
        #print(names)
        if node.dist < 0:
            node.dist = 0

        for n in names:
            if n == "Reference":
                continue
            if n not in memberships:
                memberships[n] = list()
            memberships[n].append(group)

        level_rank += 1
        group += 1

    n_ranks = list()
    for n in memberships:
        n_ranks.append(len(memberships[n]))

    max_ranks = max(n_ranks)

    for i in range(0, max_ranks):
        group = 1
        lookup = dict()
        for n in memberships:
            length = len(memberships[n])
            if i < length:
                group_id = memberships[n][i]
                if group_id not in lookup:
                    lookup[group_id] = group
                    group += 1
                memberships[n][i] = lookup[group_id]

    for n in memberships:
        count_levels = len(memberships[n])
        if count_levels < max_ranks:
            for i in range(count_levels, max_ranks):
                memberships[n].append(0)

    for n in memberships:
        print("{}\t{}".format(n, "\t".join(str(v) for v in memberships[n])))

    return memberships




def main():
    args = parse_args()
    nwk_treefile = args.in_nwk
    reference_fasta = args.reference
    vcf_file = args.in_vcf
    min_snps = args.min_snps
    min_group_size = args.min_members

    tile_size = 32
    tile_size = tile_size - 1
    mid_point = int(tile_size / 2)

    memberships = get_tree_groups(parse_tree(nwk_treefile))

    for seq_record in SeqIO.parse(reference_fasta, "fasta"):
        ref_seq = "{}".format(seq_record.seq)
        ref_length = len(ref_seq)


    reader = vcfpy.Reader.from_path(vcf_file)

    # Harcoded header
    header = ['#CHROM', 'POS', 'REF', 'ALT'] + reader.header.samples.names
    samples = reader.header.samples.names
    snp_counts = list()

    scheme = dict()

    for record in reader:
        if not record.is_snv():
            continue
        line = [record.CHROM, record.POS, record.REF]
        line += [alt.value for alt in record.ALT]

        if len(record.ALT) > 1:
            continue
        alt_base = record.ALT[0].value
        ref_base = record.REF
        position = record.POS
        count_ref = 0
        count_alt = 0
        ref_group = list()
        alt_group = list()
        tracker = 0
        for call in record.calls:
            state = int(call.data.get('GT'))
            sample_name = samples[tracker]
            if state == 0:
                count_ref+=1
                ref_group.append(sample_name)
            else:
                count_alt+=1
                alt_group.append(sample_name)
            tracker+=1

        valid_ranks = list()
        valid_groups = list()
        conflict_positions = dict()

        if(count_ref >= 1 and count_alt >= 1):
            #print(ref_group)
            ref_candidate_partitions = check_snp(ref_group,memberships,position,tile_size,ref_seq,valid_groups,valid_ranks,mid_point,ref_base,alt_base,count_ref,count_alt)
            #print(alt_group)
            alt_candidate_partitions = check_snp(alt_group,memberships,position,tile_size,ref_seq,valid_groups,valid_ranks,mid_point,ref_base,alt_base,count_ref,count_alt)
            positions = dict()



            for rank in ref_candidate_partitions:
                if not rank in scheme:
                    scheme[rank] = dict()
                for item in ref_candidate_partitions[rank]:

                    #filter out groups which are too small
                    if item['ref_count'] < min_group_size:
                        continue
                    if not item["group_id"] in scheme[rank]:
                        scheme[rank][item["group_id"]] = []
                    scheme[rank][item["group_id"]].append(item)
                    pos = item["position"]
                    if not pos in positions:
                        positions[pos] = 0
                    positions[pos] += 1


            for rank in alt_candidate_partitions:
                if not rank in scheme:
                    scheme[rank] = dict()
                for item in alt_candidate_partitions[rank]:
                    #filter out groups which are too small
                    if item['alt_count'] < min_group_size:
                        continue

                    if not item["group_id"] in scheme[rank]:
                        scheme[rank][item["group_id"]] = []
                    scheme[rank][item["group_id"]].append(item)
                    pos = item["position"]
                    if not pos in positions:
                        positions[pos] = 0
                    positions[pos] += 1
                    if positions[pos] > 1:
                        conflict_positions[pos] = item


            for rank in ref_candidate_partitions:
                for item in ref_candidate_partitions[rank]:
                    pos = item["position"]
                    #filter out groups which are too small
                    if item['ref_count'] < min_group_size:
                        continue
                    if positions[pos] == 1:
                        print("{}\t{}".format(item, "REF-Single"))
                    else:
                        print("{}\t{}".format(item, "REF-Multi"))


            for rank in alt_candidate_partitions:
                for item in alt_candidate_partitions[rank]:
                    pos = item["position"]
                    #filter out groups which are too small
                    if item['alt_count'] < min_group_size:
                        continue
                    if positions[pos] == 1:
                        #print("Success")
                        print("{}\t{}".format(item, "ALT-Single"))
                    else:
                        print("{}\t{}".format(item, "ALT-Multi"))

    #filter out groups with less than minimum number of supporting snps
    for rank in scheme:
        valid = dict()
        for group_id in scheme[rank]:
            if len(scheme[rank][group_id]) >= min_snps:
                valid[group_id] = scheme[rank][group_id]
        scheme[rank] = valid


    for n in memberships:
        heirarchy = memberships[n]
        for i in range(0,len(heirarchy)):
            if i not in scheme :
                heirarchy[i] = 0
            elif not heirarchy[i] in scheme[i]:
                heirarchy[i] = 0

        print("{}\t{}".format(n, "\t".join(str(v) for v in memberships[n])))










# call main function
if __name__ == '__main__':
    main()