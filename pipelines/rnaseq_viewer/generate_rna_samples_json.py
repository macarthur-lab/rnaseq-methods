import argparse
import collections
import logging
import os
import pprint

from sample_metadata.utils import get_joined_metadata_df

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s')
logger = logging.getLogger()
logger.setLevel(logging.INFO)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-b", "--batch-name", help="optional batch name")
    args = p.parse_args()

    return args

def main():
    args = parse_args()
    logger.info("Args:\n" + pprint.pformat(args.__dict__))

    df = get_joined_metadata_df()

    logger.info("Done")


if __name__ == "__main__":
    #main()
    pass

#%%


df = get_joined_metadata_df()

df.columns

#%%
batches = collections.Counter(df.star_pipeline_batch)

print(batches)


#%%

batch = "batch_1_muntoni"
if batch:
    df = df[df.star_pipeline_batch == batch]

#%%

df.star_pipeline_batch

df.sample_id
df.star_bam
df.star_bai

df.junctions_bed
df.coverage_bigwig
