# Changelog

## [1.0.0](https://www.github.com/clinical-genomics-uppsala/TwiMM/compare/v0.1.1...v1.0.0) (2025-11-07)


### ⚠ BREAKING CHANGES

* new input functions from HG v3.3.0 are used

### Features

* add a script for extracting translocatins from sniffles2 files ([8e63cbd](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/8e63cbd49f17f6192965a9ae151c4a4d00e42caa))
* add bcftools view filtering ([d11ccf0](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/d11ccf096c2eb9c2a5806bbb67d11c67aa6e763c))
* add claires_to rules ([43f2090](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/43f2090f173481b85724d498df956bede8e65224))
* add claires-toand connect cnvkit*rules ([42de1f2](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/42de1f2e25991d2e8d14a74b940588faa338ab40))
* add cnvkit call, vcf + vep ([6a507ef](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/6a507ef8a21a4e371e3d19636b3c95dc350f0bb2))
* add cnvkit vcf as input to rule compile_xlsx_report ([d575307](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/d575307e393a25ca2d34598de14a6e0df5c32f88))
* add cnvkit_batch as a local script ([17a2615](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/17a2615fd3acdbc1a4926e26677f75797be52a74))
* add compile_xlsx_report rule ([d952468](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/d952468c9de342b6dac2a9e09b5934a1c3e3ac2c))
* add fix_af.py to fix AF field ([cc0d3c7](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/cc0d3c735d35a0788278030f9ca492ec0e0f17b9))
* add handling of missing BAF values in cnvkit vcf files ([cd16c76](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/cd16c7674e954005e4846c66f0a79c7f593d244b))
* add local variant of cnvkit_call ([4cca549](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/4cca5496ffa43b34bd8ccdbf2e9e11f141beb0e1))
* add logging (debugging) and type hints to the functions ([137a1ea](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/137a1eaa4f77e4895f8d63b0eb1a6b8c82a37df7))
* add more columns to be extracted from SV VCF ([72d1a3e](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/72d1a3e390e7ae6f5bae37854b68b5737c077aac))
* add phasing of deepsomatic vcfs ([ca71c36](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/ca71c3601d29f7c3f65100884d164eae66baf84a))
* add script to parse SNV VCF ([b0f93d4](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/b0f93d4a6c0a6a8f0eba0e2a7753a556b41df92e))
* add tab for CNVkit VCF results ([f9e368d](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/f9e368ddb920ad4788af72f1c53b8993fefecab9))
* add translocation tabs ([7ba71f6](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/7ba71f64c3a2b5e2101758ef8d42884daf5cf087))
* make IDID variants tab in xlsx report ([2c67bf6](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/2c67bf69c6540288ab3deb34931299a66d3ecbeb))
* make the output temporary ([35cd56b](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/35cd56bc689ee445a3de75f4b8995936be8f1e04))
* make xlsx report with two tabs ([74dc66a](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/74dc66acd6cb72ce6f96799271b1874a402b8e5b))
* Merge pull request [#18](https://www.github.com/clinical-genomics-uppsala/TwiMM/issues/18) from clinical-genomics-uppsala/develop ([a614bd4](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/a614bd47e44cffb1e30736f2da71708938b69209))
* remove hificnv ([31ca1c2](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/31ca1c27ee185b7884c840f30bf1a706440f7e99))
* update cnvkit version in cnvkit call ([8f50721](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/8f507210cb2a8f88a303b849cadc1d2bbd6856f7))
* whathap_haplotag works again + filtering of sniffles results ([e714f9c](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/e714f9cd97d436a87d72f0d03ada5fa3555596b0))


### Bug Fixes

* add entries required by cnvkit/reports rules ([f181fd9](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/f181fd94969869ec0cac8e24a02e916571d89169))
* add fix_af ([aabbff4](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/aabbff4161382ed5d4c63f0016b9b5ffa1ce9f6d))
* add functioning logging ([e65b7d0](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/e65b7d06ad0722ac7ac85c9c03843eb909d2ecc0))
* add imports of the new input functions to common.smk ([b666624](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/b6666243cc1f90fa1deea2279246e2dfd47425e0))
* cnvkit batch wrapper bug ([216715e](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/216715e2404e1857a5e6a1958a527e93cc0fa5ec))
* correct clairs yaml entries ([6642fdb](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/6642fdb7fec0ac74e12d51ed5decd5b3e9ad5b03))
* correct file paths to the output of whatshap rules ([69671eb](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/69671ebf587d426e5680effe5eeafb7f92526ad3))
* correct output filenames; comment out not needed names ([74633eb](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/74633ebbb4b39d402f78d1424332c26680ff1147))
* fix  merge_cnv_json and fix_af rules ([b20c7ac](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/b20c7ac9a2324c256c97020f0c4cf5fc85e0e04d))
* hard-code inputs in cnvkit_batch ([bd76af9](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/bd76af90df2eeb75c03a03c53e68a16923ed638d))
* update Snakefile ([33eebf6](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/33eebf69498991870b73b910f1d622f42e3ee926))
* use rule fix_af from HG not the local variant ([0bb6ab5](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/0bb6ab50dcc3084bdc8dcc26fd04db8ea572bdad))


### Documentation

* add rulegraph.png with claires and correctly connected cnvkit/report rules ([12940ce](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/12940cebeb3ea1166f9ddee6afa92a6256acc0aa))
* update rulegraph ([7a892fc](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/7a892fcec38ff36d958eb11a63531200dbc42e84))
* update rulegraph ([5e80113](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/5e80113b42e55f530d6fa00469309e70cfb1d4d8))
* update rulegraph ([f998472](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/f998472ac49b4892ec2a04b98119a99b3ce3a395))
* update rulegraph ([6f2f8a6](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/6f2f8a6bd3ab514c1ef6a07079b7f18bba5d1f80))
* update rulegraph.svg ([3fdb211](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/3fdb2116befc2cafa805b7175666c21f88427be6))
* update rulegraph.svg ([8c1a08b](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/8c1a08bad7dce4f8a7f1d6f5add49dfc12656e56))
* update rulegraph.svg ([9e25595](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/9e25595ae3a90408e719493bc7cfff50e88fe04b))
* update rulegraph.svg ([b2576df](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/b2576dfc2f75d9472344c264132224ee7352d696))

### [0.1.1](https://www.github.com/clinical-genomics-uppsala/TwiMM/compare/v0.1.0...v0.1.1) (2025-08-22)


### Documentation

* add credits.md ([6ecb0de](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/6ecb0de99cf7b06c9272e50a1180a6edbb05aabe))
* add overview.md (incl. to mkdocs.yaml) ([54157ff](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/54157ffb27c5754efe4cd61445d7203b9988819f))
* add result_files.md ([bc080b4](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/bc080b491332d0fd1571450803b61e0e913c17e2))
* fix hierarchy in overview.md ([301a28a](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/301a28ad33859666f87314557e5b69ecaf7ee825))
* fix links/badges in README.md ([13204e3](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/13204e33bb0b4e3780d1eb850c1fa570572eb4cc))
* full structure of mkdocs.yaml ([b7a982f](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/b7a982fb870e89a11bd29ceccae688a5b27a083a))
* update badge groups in README.md ([f732cfc](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/f732cfc4a63b647c1a228fbc1ff01bbf580891b0))
* update docs ([596124f](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/596124f2d65fdc74747f703e079fcb9d553f0eab))
* update index.md ([a9d1fb7](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/a9d1fb7cdce0edc1c8335a9729670dea07ab0dbf))
* update index.md and intro.md ([6105415](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/6105415b25f1c21c1ae97e545e66daa49029017f))
* update intro.md and index.md ([81ab675](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/81ab675d2913a35a54d1ca1a25f2db8a7c550277))
* update README.md ([bc6db96](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/bc6db9608f88540a73b9b95b3756da691f77b0c9))
* update README.md ([5d6236b](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/5d6236b2d204f0199e8fa40462ef8f782df49969))
* update software versions and descriptions ([f94caf6](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/f94caf66c682f49ef276a93db05db5fdea00759f))
* update softwares.md ([1f58c0c](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/1f58c0c786e2381c67c73735747599b1438cca56))

## 0.1.0 (2025-04-11)


### Features

* add gzip+tabix and filtering of VCFs ([b541c27](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/b541c2782d94ea37fbb71fff74347d9383a1bae8))
* add new modules and small fixes ([946fac8](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/946fac8d74060527d6d712b7a166367b210e46c5))
* add output files entry for pbsv ([1b3fc56](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/1b3fc5653904f1349da26d2ce6fc5b8c260c9469))
* add pbsv entry ([89be387](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/89be3870b4886148d7d14f9a44d00b930be1d9c1))
* add pbsv_filtering ([ffb61bc](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/ffb61bcc1bd674ece0149732bac0113c7ccba4ca))
* add rules cnvkit_batch and cnvkit_diagram ([d05c0b4](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/d05c0b45bb815f04d6d213d9239b370a17476f3f))
* add sniffles2, remove pbsv filtering, move hificnv together with the other ruls from cnv_sv ([6e8cd32](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/6e8cd32e0c5eab46c228b4368f0dc0f86092bbd2))
* create bcftools.smk ([53fd258](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/53fd258af5c0b7923dcd4c5f263f5b8a04123ea0))
* include rule deepsomatic ([14ec726](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/14ec726176a9828e0aa331130381e968cc010b06))
* include rules/pbsv.smk ([b67ca45](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/b67ca45ba27f3bbb82bfd0efb38280da5a04e5b6))
* update rulegraph ([9f82f59](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/9f82f594dde5ed39460e407e1342702bec83a3b7))


### Bug Fixes

* add correct entry (unit_bam) to config.yaml ([6f3dcf8](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/6f3dcf8937c7f5260e768fc23905671cba02c5d1))
* add platforms PACBIO & ONT ([d1db5f6](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/d1db5f62083386bb899aa003740ff04cf2f9e035))
* add T to input and output file names ([d3ad63f](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/d3ad63f75c8c0b5977cc92e84c809d708ad1eece))
* add use hificnv from cnv_sv module ([63ad9df](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/63ad9df7e392e141efb158b5ad130fb9be57de93))
* bcftools entry in accord with the rule ([8422135](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/8422135afb49980d63fcb4fad949af11ffb79d89))
* functioning and inegrated with Hydra snakefile ([591fbe3](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/591fbe35e799a2d87429060343e2b59109d05d1e))
* input files according to what is required by Snakefile ([0aa6be9](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/0aa6be989b6ea4107707bb76aaba1b5f5b976163))
* make correct input and output file paths ([f91ed24](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/f91ed241955ec3973a610fc2c2f93db3386d3f7f))
* move multiqc & add fai, trf, design_bed to reference ([1505e24](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/1505e248600a19c8d845eb89a385a780c5f97975))
* new entry names, docker image paths ([6050885](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/6050885cd5e99ee5d7fb8d6e840dd5b06e40ff96))
* remove 'sv' from output file extension ([1f5f220](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/1f5f220e3251c7e7b2cd88868bc190793e3a001d))
* two separate entries for two pbsv rules (discover & call) ([453ee5c](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/453ee5c33aab3af5353c58b2fc9ccca57c697993))
* update gydra-egentics version to 3.0.0 ([d11db2c](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/d11db2c14a9dce035f148b37de1309440e6d39f1))
* use correct versions of cnv_sv, snv_indels and deepsomatic_tn ([7f4a0bc](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/7f4a0bc3798ad0bf98409e569d8acd4c6baef329))
* use modules snv_indels, cnv_sv from hydra ([f02348c](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/f02348c282b5c538de914ff05ef4aded1965636c))
* use rule cnvkit_batch with input from the module ([0c3e1c9](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/0c3e1c9a1318aa547e67451bd5eb7e78beaae4f3))
* use two separate pbsv rules (call & discover) instead of one ([c61c603](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/c61c60307ce95f8bca955bc186e954145cb32594))


### Documentation

* remove rulegraph.png ([91c5fed](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/91c5fed849b8b1fc8c6afadfc187b29e68f55f0c))
* update rulegraph ([4cd0966](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/4cd0966efe2e6480c779d2f8c37b37734d13998f))
* update rulegraph ([a9e2d51](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/a9e2d51ce75b8fa9145261502e3d50a08ea3bc06))
* update rulegraph ([d58299a](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/d58299aec66577fc2f653326191bde4c2ec0ade6))
* update rulegraph ([88f858c](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/88f858c01916e4d352515f95ac3e99351bca98ec))
* update rulegraph.svg ([21a6e3c](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/21a6e3c2f03230e796f859a0d7dbe390627802fc))
* update rulegraph.svg ([1b3cf35](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/1b3cf35373b50b170b6acb5ee5ea398c73262524))
* use local files as starting BAM files ([b890fbc](https://www.github.com/clinical-genomics-uppsala/twist_myeloma_pipeline/commit/b890fbc426484a7a47cb8b5073181a5a5a8d6203))
