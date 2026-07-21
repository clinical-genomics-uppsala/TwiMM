# Changelog

## [1.7.0](https://github.com/clinical-genomics-uppsala/TwiMM/compare/v1.6.0...v1.7.0) (2026-07-21)


### Features

* add bgzipped and tabixed files to the output ([ef7ea21](https://github.com/clinical-genomics-uppsala/TwiMM/commit/ef7ea2189cacf952f1096d40ad4cc1cd0d6155f7))
* add caller column to both snv and tp53 tabs ([7d73abe](https://github.com/clinical-genomics-uppsala/TwiMM/commit/7d73abe613772b73c40d0d655df2d1636e3b5f32))
* add config placeholder expansion function ([2ea34b1](https://github.com/clinical-genomics-uppsala/TwiMM/commit/2ea34b19796192fd83dedbb6c73f80157da17dba))
* add configurable CNV tab filters to xlsx report ([115880c](https://github.com/clinical-genomics-uppsala/TwiMM/commit/115880c27d133715292e8b0a240d5f39ad7aa6f5))
* add configurable SVDB database files ([d2ababf](https://github.com/clinical-genomics-uppsala/TwiMM/commit/d2ababf55765722c0c8102ffa80ed7853270c8c8))
* add cytobands ([27b7456](https://github.com/clinical-genomics-uppsala/TwiMM/commit/27b74561f95f3f2fddaaf6f9812b58a02cc35563))
* add filtering of VCF files form Severus and PBSV ([40d99be](https://github.com/clinical-genomics-uppsala/TwiMM/commit/40d99bec33588da806c116b9b34155e6560431b2))
* add rule svdb_build ([58d69e2](https://github.com/clinical-genomics-uppsala/TwiMM/commit/58d69e263da8c02b819881ad22102fcc87d94cc7))
* add rules to write caller provenance in snv vcf files and concatenate them; update snv helper function ([1472360](https://github.com/clinical-genomics-uppsala/TwiMM/commit/1472360564fda3ac54cfe9b23f90c7cc55d1f3c1))
* add section in reference with SV databases ([56f5cf3](https://github.com/clinical-genomics-uppsala/TwiMM/commit/56f5cf362f5039d6b6c1cd3c8ad1363dd6225cb6))
* add severus_t_only rule ([#66](https://github.com/clinical-genomics-uppsala/TwiMM/issues/66)) ([61f3520](https://github.com/clinical-genomics-uppsala/TwiMM/commit/61f3520614dfa6b3dd4668c1c97cf9b7ba0723cc))
* add software versions tab ([#67](https://github.com/clinical-genomics-uppsala/TwiMM/issues/67)) ([fbc4550](https://github.com/clinical-genomics-uppsala/TwiMM/commit/fbc4550e3406b2b8452bc29b07d2cdf4743f82d1))
* add svdb & pbsv rules ([973393a](https://github.com/clinical-genomics-uppsala/TwiMM/commit/973393a02a6d1b9d76265ac5dd8a49e82b4e945f))
* add vacmap aligner ([#65](https://github.com/clinical-genomics-uppsala/TwiMM/issues/65)) ([380b56c](https://github.com/clinical-genomics-uppsala/TwiMM/commit/380b56c375a772ccc5037489d523513e10c74ce0))
* bgzip and tabix the final svdb query vcf file ([5ed5173](https://github.com/clinical-genomics-uppsala/TwiMM/commit/5ed517352bf83747464d15163b5cfac068da40ad))
* cnvkit pon with bigger bin size (2kb) ([a2df88c](https://github.com/clinical-genomics-uppsala/TwiMM/commit/a2df88cc7d101271ac1fff3795c63b9bf598bf8b))
* hard-filter CNVkit germline VCF ([584b624](https://github.com/clinical-genomics-uppsala/TwiMM/commit/584b6245a216cd92ca348d26b83c4a319a44fecf))
* introduce updates compatible with HG v4 ([60e77cb](https://github.com/clinical-genomics-uppsala/TwiMM/commit/60e77cbecf3fdbb0c4a710ac9fc9006edb38d2fc))
* log edge case when neither VAF or AF is present ([5967aa9](https://github.com/clinical-genomics-uppsala/TwiMM/commit/5967aa96b0507b968147dea205c02fcfd5a3f62e))
* SVDB results to Excel report ([#64](https://github.com/clinical-genomics-uppsala/TwiMM/issues/64)) ([168f009](https://github.com/clinical-genomics-uppsala/TwiMM/commit/168f0099d09b771b8ae2b3ff589bb4e8a4ac2748))
* upd config schema file with sv_databases and svdb_build ([bcc6af2](https://github.com/clinical-genomics-uppsala/TwiMM/commit/bcc6af2a4a25e55f94b46da03c764039c9fb692a))
* update config.schema.yaml with new keys ([ccf7dc4](https://github.com/clinical-genomics-uppsala/TwiMM/commit/ccf7dc4225a453807a55ebe4dd59fe2fd677cc11))
* update filter entries ([2102aac](https://github.com/clinical-genomics-uppsala/TwiMM/commit/2102aac002053698791f31e2671ca49f71bf7120))
* use 2 aligners configuration ([f666370](https://github.com/clinical-genomics-uppsala/TwiMM/commit/f6663706b1e99319cb61a5be429712a62932f401))
* use hlper funciton to create db_string parameter for svdb query ([d86cb7e](https://github.com/clinical-genomics-uppsala/TwiMM/commit/d86cb7e5c1c5f628fcce933a53a8f8353b61d433))


### Bug Fixes

* add additional design bed file ([34938ea](https://github.com/clinical-genomics-uppsala/TwiMM/commit/34938ea39cad43b4e73dea72886dd8d4971b516b))
* add benchmark, threads and other properties to svdb_build entry ([eb7d344](https://github.com/clinical-genomics-uppsala/TwiMM/commit/eb7d344f33d80947040470f6fc68e34c2dd5677e))
* add configs for cnv and germline filters ([2f55222](https://github.com/clinical-genomics-uppsala/TwiMM/commit/2f552229ec6db5b9cd6001e8909931b2f42f48a6))
* add cytoband file to integration test ([7413d98](https://github.com/clinical-genomics-uppsala/TwiMM/commit/7413d98ae3c56c55c9fe1ea6921c11520eadc84c))
* add db prefix to params ([ca96fe2](https://github.com/clinical-genomics-uppsala/TwiMM/commit/ca96fe23076b36b999d302da47547fbe7d59d247))
* add more files for integration test ([2a27dc3](https://github.com/clinical-genomics-uppsala/TwiMM/commit/2a27dc3e5610d20bc2d404bccc2a656dcbc31480))
* add vcf report config to integration test ([11f13e2](https://github.com/clinical-genomics-uppsala/TwiMM/commit/11f13e293bacbf4003bff7c272933967c693c897))
* additional VAF guard for edge case ([43bba38](https://github.com/clinical-genomics-uppsala/TwiMM/commit/43bba38c13efe9b2194c8f636648d27f9d536df0))
* cap SV/Translocation ALT column width via config ([ab0dac2](https://github.com/clinical-genomics-uppsala/TwiMM/commit/ab0dac2a5a300992bbb1379d100e186c16ba4557))
* compatibility with python 3.9 ([0c5dc20](https://github.com/clinical-genomics-uppsala/TwiMM/commit/0c5dc202bb4ec684538e5396d9847bca378ef1b9))
* correct dual-caller concat pipeline bugs ([4a48213](https://github.com/clinical-genomics-uppsala/TwiMM/commit/4a482132443ea78310f3ed5fd3a0102cd5f44ecd))
* correct flag in svdb query when using with sqlite ([05691d4](https://github.com/clinical-genomics-uppsala/TwiMM/commit/05691d41c4cb0657869352dd5f4b74fea0f2ae72))
* derive PBSV VAF from AD and format all SV VAF to 2 decimal places ([8ebff25](https://github.com/clinical-genomics-uppsala/TwiMM/commit/8ebff25c172ce70f027b1fe10ce29d60b5c5c9f2))
* derive svdb_build prefix from output via lambda to fix snakefmt lint. ([aa90daa](https://github.com/clinical-genomics-uppsala/TwiMM/commit/aa90daa7add0da6d6bb5e13695583d0a50e6dddc))
* generate real intervals list for picard ([cebe544](https://github.com/clinical-genomics-uppsala/TwiMM/commit/cebe54470037ee70345cfe960ce1b871a4ed7a10))
* handle DeepSomatic VAF field in fix_af and parse_vcf_line ([21c44b3](https://github.com/clinical-genomics-uppsala/TwiMM/commit/21c44b368122bb4d760f0f0b8592732100eceab2))
* make hydra-genetics v4 default environment ([d762d9e](https://github.com/clinical-genomics-uppsala/TwiMM/commit/d762d9e137352c56fcea90c5ad59699837aada11))
* move CALLER to column in SNV/TP53 tabs and preserve headers on empty tabs ([952eff4](https://github.com/clinical-genomics-uppsala/TwiMM/commit/952eff42724cd2f023349549152133ce369d1ac5))
* normalize caller to a consistent case ([fd51f1e](https://github.com/clinical-genomics-uppsala/TwiMM/commit/fd51f1eeaab8e25c8334eb03db6bea25fce7c919))
* populate SUPPORT from variant read depth for Severus and PBSV records ([c02b202](https://github.com/clinical-genomics-uppsala/TwiMM/commit/c02b2023884f3023e05061ea7c3cb808181e0818))
* rename path to VEP's cache dummy file ([ee487cc](https://github.com/clinical-genomics-uppsala/TwiMM/commit/ee487cc0cd672a2e10a01a0c3c0eb686c9bd5a66))
* safe convert coverage to int ([5f6e00f](https://github.com/clinical-genomics-uppsala/TwiMM/commit/5f6e00f7c5d33fd161f48562dcda7974cb6814e4))
* sanitize slashes in pipeline version string before use in file paths ([d35cf8d](https://github.com/clinical-genomics-uppsala/TwiMM/commit/d35cf8dfab1eba26218c84d45d33aed9467eb512))
* set wildcard constraints to T only ([4facc1a](https://github.com/clinical-genomics-uppsala/TwiMM/commit/4facc1ab1e307b62613009889db2a87d2baea13e))
* trim cytobands file ([5c9e1cb](https://github.com/clinical-genomics-uppsala/TwiMM/commit/5c9e1cbf03c312f5bcaf66327f344e4cbbae2098))
* unwrap pysam SVLEN tuple to prevent chr14 SVs being dropped ([51c4cd5](https://github.com/clinical-genomics-uppsala/TwiMM/commit/51c4cd53648052ca7782055282211e6516981cc6))
* update config with all the required paths and placeholders ([4d274c7](https://github.com/clinical-genomics-uppsala/TwiMM/commit/4d274c7d6a06a444a11fb46ef393a13b1d759932))
* update lock file ([1ce4318](https://github.com/clinical-genomics-uppsala/TwiMM/commit/1ce4318b2fa6f9598d9c751e6151fc6547b56274))
* update path to correct gnomad sv file ([407116e](https://github.com/clinical-genomics-uppsala/TwiMM/commit/407116ea8596cb4479013d6d666e3e1b275c594a))
* update paths to ref files ([039fb84](https://github.com/clinical-genomics-uppsala/TwiMM/commit/039fb84449edea790681d22d6ea05246db9ec491))
* update the handling around raw_cov to not iterate strings/Scalars ([3105016](https://github.com/clinical-genomics-uppsala/TwiMM/commit/31050165f203a15c00f144a122b46efb4b6ac44f))
* use a helper function instead of hard-coded input paths ([422ab06](https://github.com/clinical-genomics-uppsala/TwiMM/commit/422ab06ded6eb3ce7f60a04de930ae9216346d55))
* use direct key access in svdb_build input to catch missing key ([3a9455c](https://github.com/clinical-genomics-uppsala/TwiMM/commit/3a9455cabc9ec144adbff6d8a90a1b739cef5300))
* **xlsx-report:** fix SV VCF parsing and add population frequency columns ([2800366](https://github.com/clinical-genomics-uppsala/TwiMM/commit/28003660047ace656d3a307dabfd781c458ac481))
* **xlsx-report:** harden pysam INFO field handling and remove dead code ([8b9f711](https://github.com/clinical-genomics-uppsala/TwiMM/commit/8b9f7110ddc567117b7db67b5e534b1c2b9c198e))


### Reverts

* remove svdb build ([4cbf02a](https://github.com/clinical-genomics-uppsala/TwiMM/commit/4cbf02a3782a049f344ab573a49c10167c975f5f))


### Documentation

* clarify commemts on filter usage ([77d4e60](https://github.com/clinical-genomics-uppsala/TwiMM/commit/77d4e60c4c857505ee7fef8410d828a80878ab24))
* fix a docstring typo ([6b96e63](https://github.com/clinical-genomics-uppsala/TwiMM/commit/6b96e637ebacbea0e2104323d984484a5a9d88d3))
* fix typos ([930e1b0](https://github.com/clinical-genomics-uppsala/TwiMM/commit/930e1b0b3a5c992c2a131ba6e419c380a71a776e))
* material theme ([3c80db3](https://github.com/clinical-genomics-uppsala/TwiMM/commit/3c80db36a34cffb8d7650343f469919396e1b984))
* remove deepsomatic rulegraphs ([f76d9a5](https://github.com/clinical-genomics-uppsala/TwiMM/commit/f76d9a5a160fc0757355535db9c6bf0c823228a4))
* revise and update output files list ([b80010c](https://github.com/clinical-genomics-uppsala/TwiMM/commit/b80010ce54c020928e6d40e95e461869cfee9697))
* unify output file paths ([b72538b](https://github.com/clinical-genomics-uppsala/TwiMM/commit/b72538b641aeb85614d6efbaca46b78443685b87))
* upd hydra and snakemake version in softwares.md ([1032cf6](https://github.com/clinical-genomics-uppsala/TwiMM/commit/1032cf68068cd8eb2dea73d08f56ff13066afdd0))
* upd overview.md ([e80863e](https://github.com/clinical-genomics-uppsala/TwiMM/commit/e80863ea4f1bc72b7541d6079344ce30cd09f652))
* upd README.md ([5e6e1ac](https://github.com/clinical-genomics-uppsala/TwiMM/commit/5e6e1ac71c61ab0120245028392562f5cfdc8e41))
* upd rulegraph ([6768f5f](https://github.com/clinical-genomics-uppsala/TwiMM/commit/6768f5f47d9672e3b0a8363ac0234f7009ae3167))
* update documentation files ([d17fa4c](https://github.com/clinical-genomics-uppsala/TwiMM/commit/d17fa4cdbefb0c0d5bcffedae0a5a70a6c8fb991))
* update gnomad sv filename ([391b33a](https://github.com/clinical-genomics-uppsala/TwiMM/commit/391b33a79fc75c91529ff691dca8d6c90528f164))
* update README.md ([546cccf](https://github.com/clinical-genomics-uppsala/TwiMM/commit/546cccf40cead26ffe94d4bb12473664dd6332cb))
* update rulegraph/dot files ([a27dd96](https://github.com/clinical-genomics-uppsala/TwiMM/commit/a27dd968d8e3123ac76f5162958beac303fd73cf))
* update rulegraphs ([83868c0](https://github.com/clinical-genomics-uppsala/TwiMM/commit/83868c0b203b207e586f48700063e995510ee7e0))
* update rulegraphs ([5bf865b](https://github.com/clinical-genomics-uppsala/TwiMM/commit/5bf865bd562c5df2979cc50d8a77cac9bcc4db3d))

## [1.5.0](https://github.com/clinical-genomics-uppsala/TwiMM/compare/v1.4.0...v1.5.0) (2026-02-11)


### Features

* add function to dyna,ically append ref genome to VEP extra in the config; fix datestring as pipeline version ([07a9393](https://github.com/clinical-genomics-uppsala/TwiMM/commit/07a9393a425d67abb7f0a4f50747f8690fab2978))
* add picard rules ([7d739be](https://github.com/clinical-genomics-uppsala/TwiMM/commit/7d739be89e7a4df5f8b468bd3579a964725544dc))
* Introduce configurable pipeline version and dynamic VEP FASTA path ([dcf12fd](https://github.com/clinical-genomics-uppsala/TwiMM/commit/dcf12fd9512b54e6dc8a29b8f521976209570a41))


### Bug Fixes

* pin reports version to a commit before v0.11.0 ([63f7a1a](https://github.com/clinical-genomics-uppsala/TwiMM/commit/63f7a1a9ad425c5ed2daa64c70ba1940444922aa))
* use ref_genes in merge_cnv_json ([00e9004](https://github.com/clinical-genomics-uppsala/TwiMM/commit/00e9004c482a0bc1d55fd533a0cc7cb53689e6f7))


### Reverts

* remove ref_genes from merge_cnv_json config and schema ([47075c1](https://github.com/clinical-genomics-uppsala/TwiMM/commit/47075c18ad51de9988777f4568634d81a6d4b468))


### Documentation

* rm old rulegraph ([d8440eb](https://github.com/clinical-genomics-uppsala/TwiMM/commit/d8440ebd812e77757a3918b7fa0548b348701059))
* upd rulegraphs and README.md ([e214d87](https://github.com/clinical-genomics-uppsala/TwiMM/commit/e214d87a8e3c9fb5f3a615b7c4356fdc3b89cbf1))
* update rulegraphs ([7e65f87](https://github.com/clinical-genomics-uppsala/TwiMM/commit/7e65f87057b246ffbab6e1fd941b7962aa612aae))

## [1.4.0](https://github.com/clinical-genomics-uppsala/TwiMM/compare/v1.3.0...v1.4.0) (2026-01-22)


### Features

* add filtering of genes with SNVs ([6166c06](https://github.com/clinical-genomics-uppsala/TwiMM/commit/6166c063d5d22533e6d2d79014286a4128a19188))
* add genes for filtering as param ([292603f](https://github.com/clinical-genomics-uppsala/TwiMM/commit/292603fb31dad8ffa09289fea42868e4b0a16c15))
* add matching column width ([d837221](https://github.com/clinical-genomics-uppsala/TwiMM/commit/d83722131ab05a7571b2c9a33d5edf3d7b8427e2))
* add MB to the column header in CNV tab ([ea48652](https://github.com/clinical-genomics-uppsala/TwiMM/commit/ea48652b018780e4faea2e3b282fa10b48a56777))
* use rule mosdepth instead of mosdepth_bed ([fc4c568](https://github.com/clinical-genomics-uppsala/TwiMM/commit/fc4c568be8f53091bce2d67ed7c7216e425cde73))


### Bug Fixes

* update version_files array ([fcf8d44](https://github.com/clinical-genomics-uppsala/TwiMM/commit/fcf8d449a02d44933e2845596003ee8d1e3ad527))

## [1.3.0](https://www.github.com/clinical-genomics-uppsala/TwiMM/compare/v1.2.1...v1.3.0) (2025-12-15)


### Features

* add sequali and multiqc_long_read ([62c2f14](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/62c2f14e29cbec013a75ee1e16df4c8e9d6f0cf9))

### [1.2.1](https://www.github.com/clinical-genomics-uppsala/TwiMM/compare/v1.2.0...v1.2.1) (2025-12-10)


### Documentation

* fix formatting of the pippeline overview ([deac0d7](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/deac0d752b95333cc815f0a76227c771d53e9ef0))
* update overview.md ([d13ef14](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/d13ef14224af85da0dc51eac93f17a05c17aa887))
* update pipeline description ([8af888a](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/8af888a134861287335aab70491ee5e87cdc2835))
* update site name ([daa7820](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/daa7820d2f8d37b4f05bf93d824b82500c20960b))
* update the documentation ([90b8000](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/90b8000696b3c1dc5c8eb047c75da3ec34387c7e))

## [1.2.0](https://www.github.com/clinical-genomics-uppsala/TwiMM/compare/v1.1.0...v1.2.0) (2025-12-05)


### Features

* add a function determining which snv caller to use ([e6fc7e7](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/e6fc7e720ef97375cc854540fc9188f15d053c9d))
* add conditional executin of clairs or deeps ([6d88433](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/6d88433535bb3ec6aa0cd7ce307beed6beb1df90))


### Bug Fixes

* chnage germline_vcf filename ([590fc31](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/590fc311fb229a5de0d58c80bcad764c9abcdcc1))
* remove f-strings in get_snv_caller_output ([b499374](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/b499374d99b0504f2ac8dfa68b21cf9b76af27be))
* rename SYMBOL to GENE ([3dee681](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/3dee6815ebf9b54be3d2959766f41bb3b935356f))
* return PHENO & SOMATIC to vep_fields ([7aaf6e4](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/7aaf6e448a962419a7dbe20fc4e4d1484c56ef27))


### Documentation

* add two dags ([c537dce](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/c537dce05f2ac7e95230268475e08b5acd1d9eb0))
* add two DAGs ([2712b8a](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/2712b8aafd2e13d0a0d58aa03b2bddd2d114926b))

## [1.1.0](https://www.github.com/clinical-genomics-uppsala/TwiMM/compare/v1.0.0...v1.1.0) (2025-11-17)


### Features

* add vcf fileds to use by reports script ([ba891d9](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/ba891d949bdeb53113be4758d3d7c732736c8d80))
* read constants from snakemake.params ([c455146](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/c45514675226e62d173bdde4c299146a69da4521))


### Bug Fixes

* add safe conversion of empty strings to NaN ([cba8fd8](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/cba8fd837da2a14d5f81e11739560457cfdce103))
* use snakemake.config.get not snakemake.params.get ini reports script ([a133986](https://www.github.com/clinical-genomics-uppsala/TwiMM/commit/a13398699568b713b4bec61aa39da7ee5c85a5d0))

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
