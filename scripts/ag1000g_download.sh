#!/bin/bash
agpath=vo_agam_release/v3
mkdir -p $agpath/metadata
gsutil -m rsync -r gs://vo_agam_release/v3/metadata/ $agpath/metadata
mkdir -p $agpath/snp_genotypes/all
gsutil -m rsync -r -x '.*/calldata/(AD|GQ|MQ)/.*' gs://vo_agam_release/v3/snp_genotypes/all $agpath/snp_genotypes/all
echo 'ok' > 'download_complete.txt'
