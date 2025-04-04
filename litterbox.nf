#! /usr/bin/env nextflow
nextflow.enable.dsl=2

if (!params.output_directory){
    error("Output directory is required.")
}
if (!params.cat_gff3){
    error("cat_gff3 is required")
}
if (!params.cat_fasta){
    error("cat_fasta is required")
}
if (!params.hg38_fasta){
    error("hg38_fasta is required")
}
if (params.scavenge_ens){
    if (!params.ens_gtf){
        error("ens_gtf is required if scavenging an Ensembl annotation")
    }
    if (!params.ens_id){
        error("ens_id is required if scavenging an Ensembl annotation")
    }
    if (!params.chain){
        if (!params.ens_fasta){
            error("ens_fasta required if using the same assembly version Ensembl annotation")
        }
    }
    else{
        if (!params.ucsc_fasta){
            error("ucsc_fasta required if lifting from a different Ensembl version annotation")
        }
    }
    
}

process get_allowlist{
    input:
    tuple path(gencode), path(hgnc_ens)
    
    output:
    tuple path("gencode.genes"), 
    path("gencode.mito"),
    path("gencode.mito.gtf"),
    path("gencode.tx"),
    path("gencode.tx2gene")
    
    script:
    def chrM=params.mito_name_hg38
    """
    ${baseDir}/scripts/get_allowlist.py --gtf ${gencode} \
        -he ${hgnc_ens} -o gencode --mito ${chrM}
    """
}

process filter_human_annotation{
    input:
    tuple path(human_annotation),
    path(genes),
    path(tx),
    path(tx2gene),
    path(mito)
   
    publishDir "${params.output_directory}", mode: "copy"
     
    output:
    path("human.gtf")
    
    script:
    """
    ${baseDir}/scripts/filter_human_annotation.py --gtf ${human_annotation} \
        --allowed_genes ${genes} \
        --allowed_tx ${tx} \
        --tx2gene ${tx2gene} \
        --mito ${mito} \
        -o human.unsorted.gtf
    ${baseDir}/scripts/sort_annotation.py -i human.unsorted.gtf -o human.gtf
    """
}

process cleanup_cat_annotation{
    
    input:
    tuple path(cat_gff3),
    path(genes),
    path(mito),
    path(mito_gtf),
    path(tx),
    path(tx2gene)
    
    output:
    tuple path("*.main.gtf"),
    path("*.rm.genes"),
    path("*.rm.genes.notx"),
    path("*.rm.tx")
     
    script:
    def catBase = ( cat_gff3.toString() =~ /(.+)\.gff3(\.gz)?/ )[0][1]
    def scriptExtra=""
    if (params.excl_contigs){
        scriptExtra += " --excl_contigs " + params.excl_contigs
    }
    if (params.rename_contigs){
        scriptExtra += " --rename_contigs " + params.rename_contigs
    }
    if (params.drop_denovo){
        scriptExtra += " --drop_denovo"
    }
    """
    ${baseDir}/scripts/cleanup_cat_annotation.py --gff3 ${cat_gff3} \
        -a gencode -o ${catBase}${scriptExtra}
    """
}

process mito_liftoff{
    
    input:
    tuple path(cat_fasta),
    path(hg38_fasta),
    path(genes),
    path(mito),
    path(mito_gtf),
    path(tx),
    path(tx2gene)
    
    output:
    path("gencode.mito.liftoff.gtf")
    
    script:
    def cat_mito = params.cat_mito
    """
    ${baseDir}/scripts/mito_liftoff.py \
        --cat_fasta ${cat_fasta} \
        --hg38_fasta ${hg38_fasta} \
        --cat_mito ${cat_mito} \
        --allowlist_base gencode \
        -o gencode
    """
    
}

process get_ens_ids{
    input:
    val(species_id)
    
    output:
    path("*_ids.txt")
    
    script:
    """
    ${baseDir}/scripts/get_ensembl_names.R ${species_id} > ${species_id}_ids.txt
    """
}

process get_homology{
    input:
    tuple val(species_id),
    path(ids),
    path(notx)
    
    output:
    file("*.homology")

    script:
    """
    ${baseDir}/scripts/get_homologous_genelist.R ${ids} ${notx} > ${species_id}.homology
    """
}

process map_between_assemblies{
    input:
    tuple path(ens_fasta), path(ucsc_fasta)
    
    output:
    path("gencode.assmap")

    script:
    """
    samtools faidx ${ens_fasta}
    samtools faidx ${ucsc_fasta}
    ${baseDir}/scripts/map_between_assemblies.py ${ens_fasta}.fai ${ucsc_fasta}.fai > gencode.assmap
    """
}

process scavenge_ensembl{
    input:
    tuple path(homology),
    path(ens_gtf),
    path(hgnc_ens),
    path(ensmap)

    output:
    tuple path("gencode.rescued.gid"),
    path("gencode.rescued.gtf")
    
    script:
    """
    ${baseDir}/scripts/scavenge_old_ensembl_annotation.py \
        -G ${homology} \
        -g ${ens_gtf} \
        --hgnc_ens ${hgnc_ens} \
        --idmap ${ensmap} \
        -o gencode
    """
}

process lift_ensembl{
    input:
    tuple path(homology),
    path(ens_gtf),
    path(chain),
    path(hgnc_ens),
    path(assmap)
        
    output:
    tuple path("gencode.rescued.gid"),
    path("gencode.rescued.gtf")
     
    script:
    """
    ${baseDir}/scripts/lift_old_ensembl_annotation.py -G ${homology} \
        -g ${ens_gtf} \
        -c ${chain} \
        --hgnc_ens ${hgnc_ens} \
        --idmap ${assmap} \
        -o gencode
    """
}

process sort_annotation{
    input:
    tuple path(main_gtf),
    path(mito_gtf)

    publishDir "${params.output_directory}", mode: "copy"
 
    output:
    path("litterbox.gtf")
  
    script:
    """
    ${baseDir}/scripts/sort_annotation.py -o litterbox.gtf \
        -i ${main_gtf} ${mito_gtf}
    """
}

process sort_annotation_rescued{
    input:
    tuple path(main_gtf),
    path(mito_gtf),
    path(rescue_gtf)
    
    publishDir "${params.output_directory}", mode: "copy"
    
    output:
    path("litterbox.gtf")
    
    script:
    """
    ${baseDir}/scripts/sort_annotation.py -o litterbox.gtf \
        -i ${main_gtf} ${mito_gtf} ${rescue_gtf}
    """
}

process publish_norescue{
    input:
    tuple path(rmgenes),
    path(rmnotx),
    path(rmtx)
    
    publishDir "${params.output_directory}", mode: "copy"
    
    output:
    path("litterbox.genes.rm")
    
    script:
    """
    cp ${rmgenes} litterbox.genes.rm
    """
}

process publish_rescue{
    input:
    tuple path(rmgenes),
    path(rmnotx),
    path(rmtx),
    path(rescue_gid),
    path(rescue_gtf)
    
    publishDir "${params.output_directory}", mode: "copy"
    
    output:
    path("litterbox.genes.rm")
     
    script:
    """
    cat ${rmgenes} | grep -v -f ${rescue_gid} > litterbox.genes.rm
    """
} 

process combine_annotations{
    input:
    tuple path(gtf1),
    path(gtf2)
    
    output:
    path("combined.gtf")
    
    script:
    """
    cat ${gtf1} ${gtf2} > combined.gtf
    """
}    

process filter_rm{
    input:
    tuple path(rm),
    path(rescue),
    path(rescuegtf)

    output:
    path("genes.rm")
    
    script:
    """
    cat ${rm} | grep -v -f ${rescue} > genes.rm
    """
}

process rename_rm{
    input:
    path(rm)
    
    output:
    path("genes.rm")
    
    script:
    """
    cp ${rm} genes.rm
    """
}

process cat_gtf{
    input:
    tuple path(gtf1), path(gtf2)
    
    output:
    path("concat.gtf")
    
    script:
    """
    cat ${gtf1} ${gtf2} > concat.gtf
    """
}

process cat_gtf_rescue{
    input:
    tuple path(gtf1), path(gtf2), path(gtf3)
    
    output:
    path("concat.gtf")
    
    script:
    """
    cat ${gtf1} ${gtf2} ${gtf3} > concat.gtf
    """
}

process svgenes{
    input:
    tuple val(gtf1name),
    val(gtf1sub),
    path(gtf1),
    val(gtf2name),
    val(gtf2sub),
    path(gtf2),
    path(fa1),
    path(fa2)
    
    output:
    tuple path("synteny1.pdf"),
    path("synteny1.${gtf1name}.gz"),
    path("synteny1.${gtf2name}.gz"),
    path("synteny1.removed.${gtf1sub}.genes"),
    path("synteny1.removed.${gtf2sub}.genes"),
    path("synteny1.remap.segment.*.bed"),
    path("synteny1.remap.segment.*.gtf")

    script:
    def s1 = (fa1.toString() =~ /([^\/]+).(fa|fasta)(.gz)?/)[0][1]
    def s2 = (fa2.toString() =~ /([^\/]+).(fa|fasta)(.gz)?/)[0][1]
    
    """
    samtools faidx ${fa1}
    samtools faidx ${fa2}
    ${baseDir}/svgenes/svgenes.py \
        --gene_id2 ${params.gene_id} \
        --gene_name2 ${params.gene_name} \
        --species1 ${s1} \
        --species2 ${s2} \
        --fai1 ${fa1}.fai \
        --fai2 ${fa2}.fai \
        -o synteny1 \
        -1 ${gtf1} \
        -2 ${gtf2} \
        -F 
    """
}

process svgenes2{
    input:
    tuple val(gtf1name),
    val(gtf1sub),
    path(gtf1),
    val(gtf2name),
    val(gtf2sub),
    path(gtf2),
    path(fa1),
    path(fa2)
    
    output:
    tuple path("synteny.pdf"),
    path("synteny.${gtf1name}.gz"),
    path("synteny.${gtf2name}.gz"),
    path("synteny.removed.${gtf1sub}.genes"),
    path("synteny.removed.${gtf2sub}.genes")

    script:
    def s1 = (fa1.toString() =~ /([^\/]+).(fa|fasta)(.gz)?/)[0][1]
    def s2 = (fa2.toString() =~ /([^\/]+).(fa|fasta)(.gz)?/)[0][1]
    
    """
    samtools faidx ${fa1}
    samtools faidx ${fa2}
    ${baseDir}/svgenes/svgenes.py \
        --gene_id2 ${params.gene_id} \
        --gene_name2 ${params.gene_name} \
        --species1 ${s1} \
        --species2 ${s2} \
        --fai1 ${fa1}.fai \
        --fai2 ${fa2}.fai \
        -o synteny \
        -1 ${gtf1} \
        -2 ${gtf2}
    """
}

process unzip_or_passthrough_fa{
    input:
    path(fasta)
    
    output:
    path("*unzip.fa")
    
    script:
    def fabase = (fasta.toString() =~ /(.+).(fa|fasta)(.gz)/)[0][1]
    """
    if [ \$( echo -e "${fasta}" | grep ".gz" | wc -l ) -gt 0 ]; then
        gunzip -c ${fasta} > ${fabase}.unzip.fa
    else
        cp ${fasta} ${fabase}.unzip.fa
    fi
    """
}

process remap_segment{
    input:
    tuple val(id),
    path(bed),
    path(gtf),
    path(fasta)
    
    output:
    tuple path("remap.${id}.gtf"),
    path("remap.${id}.genes")
    
    script:
    """
    bedtools getfasta -fi ${fasta} -bed ${bed} -fo remap.fa
    liftoff -g ${gtf} -o remap.${id}.unfilt.gtf -u unmapped.txt remap.fa ${fasta}
    start_pos=\$( cat ${bed} | tail -1 | cut -f2 )
    chrom=\$( cat ${bed} | tail -1 | cut -f1 )
    ${baseDir}/scripts/filter_liftoff_gtf.py -g remap.${id}.unfilt.gtf -c \$chrom \
        -p \$start_pos --gene_name ${params.gene_name} --gene_id ${params.gene_id} \
        -o remap.${id}.gtf -G remap.${id}.genes
    """
}

process cat_remap_gtfs{
    input:
    path(gtfs)
 
    output:
    path("remap.all.gtf")
    
    script:
    """
    cat ${gtfs} > remap.all.gtf
    """
}

process cat_remap_gtfs2{
    input:
    tuple path(remap_cat), path(svgenes_gtf)
    
    output:
    path("remap.concat.gtf")
    
    script:
    """
    cat ${remap_cat} > remap.concat.gtf
    zcat ${svgenes_gtf} >> remap.concat.gtf
    """
}

process cat_remap_genes{
    input:
    path(genes)
    
    output:
    path("remap.all.genes")
    
    script:
    """
    cat ${genes} > remap.all.genes
    """
}

process finalize{
    input:
    tuple path(rm_orig),
    path(svgenes_pdf),
    path(rm_svgenes1),
    path(svgenes_filtered),
    path(rm_svgenes),
    path(rescue_genes)
    
    publishDir "${params.output_directory}", mode: 'copy'
    
    output:
    tuple path("litterbox.gtf"), path("litterbox.genes_removed.txt"), path("synteny.pdf")

    script:
    """
    cp ${svgenes_pdf} synteny.pdf
    cat ${rm_orig} | cut -f2 > rmnames
    wc -l rmnames
    wc -l ${rescue_genes} 
    wc -l ${rm_svgenes1}
    wc -l ${rm_svgenes}
    cat rmnames ${rm_svgenes1} | grep -v -f ${rescue_genes} | sort | uniq > rm1
    cat rm1 ${rm_svgenes} | sort | uniq > litterbox.genes_removed.txt
    ${baseDir}/scripts/sort_annotation.py -o litterbox.gtf \
        -i ${svgenes_filtered}
    """ 
}

workflow{
    
    gencode_dat = Channel.fromPath(params.gencode).combine(Channel.fromPath(params.hgnc_ens)) | get_allowlist
    if (params.human_annotation){
        Channel.fromPath(params.human_annotation).combine(gencode_dat).map{ tup -> 
            return [tup[0], tup[1], tup[4], tup[5], tup[2]]
        } | filter_human_annotation
    } 
    cat_dat = Channel.fromPath(params.cat_gff3).combine(gencode_dat) | cleanup_cat_annotation    
    mito_gtf = Channel.fromPath(params.cat_fasta).combine(Channel.fromPath(params.hg38_fasta))\
        .combine(gencode_dat) | mito_liftoff
    
    to_sort = cat_dat.map{ tup ->
        tup[0] 
    }.combine(mito_gtf)
        
    if (params.scavenge_ens){
        def fn = "${baseDir}/data/${params.ens_id}_ids.txt"
        if (!file(fn).exists()){
            ids = Channel.value(params.ens_id) | get_ens_ids
        }
        else{
            ids = Channel.fromPath("${baseDir}/data/${params.ens_id}_ids.txt")
        }
        homology = ids.combine(cat_dat).map{ tup -> 
            [params.ens_id, tup[0], tup[3]]
        } | get_homology
        
        if (!params.chain){
            // Same genome assembly
            assmap = Channel.fromPath(params.ens_fasta).combine(Channel.fromPath(params.cat_fasta)) | map_between_assemblies 
            rescued = homology.combine(Channel.fromPath(params.ens_gtf)).combine(Channel.fromPath(params.hgnc_ens)).combine(assmap) | scavenge_ensembl
        }
        else{
            assmap = Channel.fromPath(params.ucsc_fasta).combine(Channel.fromPath(params.cat_fasta)) | map_between_assemblies 
            rescued = homology.combine(Channel.fromPath(params.ens_gtf)).combine(Channel.fromPath(params.chain)).combine(Channel.fromPath(params.hgnc_ens)).combine(assmap) | lift_ensembl
            
        }
        if (params.synteny){
            genesrm = cat_dat.map{ tup -> tup[1] }.combine(rescued) | filter_rm
            gtfcat = to_sort.combine(rescued.map{ tup -> 
                tup[1]
            }) | cat_gtf_rescue
        }
        else{
            to_sort.combine(rescued.map{ tup -> 
                tup[1] 
            }) | sort_annotation_rescued
            
            cat_dat.map{ tup -> 
                [tup[1], tup[2], tup[3]]
            }.combine(rescued) | publish_rescue        
        }
             
    }
    else{
        if (params.synteny){
            genesrm = cat_dat.map{ tup -> tup[1] } | rename_rm
            gtfcat = cat_gtf(to_sort)   
        }
        else{
            sort_annotation(to_sort)
            cat_dat.map{ tup ->
                [tup[1], tup[2], tup[3]]
            } | publish_norescue
        }
    }
    
    if (params.synteny){ 
        // Now run svgenes
        gtf1 = Channel.fromPath(params.gencode).map{ gtf -> 
            def match = ( gtf =~ /([^\/]+).(gtf|gff|gff3)/)[0]
            return [ match[0], match[1], gtf ]
        }
        gtf2 = gtfcat.map{ gtf -> 
            def match = ( gtf =~ /([^\/]+).(gtf|gff|gff3)/)[0]
            return [match[0], match[1], gtf ]
        }
        
        // Contents:
        // Species1 (human) GTF
        // Species2 (other) GTF
        // Species1 (human) removed gene list
        // Species2 (other) removed gene list
        // list of BED segments for re-mapping removed species2 genes
        // list of GTF genes for re-mapping removed species2 genes
        svgenes_out = gtf1.combine(gtf2).combine(Channel.fromPath(params.hg38_fasta)).combine(Channel.fromPath(params.cat_fasta)) | svgenes
        
        // Get all segments to re-map 
        remapped = svgenes_out.flatMap{ tup -> 
            tup[5]
        }.map{ fn -> 
            def uid = (fn =~ /synteny1.remap.segment.([0-9]+).bed/)[0][1]
            return [uid, fn]    
        }.join(svgenes_out.flatMap{ tup -> 
            tup[6]
        }.map{ fn -> 
            def uid = (fn =~ /synteny1.remap.segment.([0-9]+).gtf/)[0][1]
            return [uid, fn]
        }).combine(unzip_or_passthrough_fa(Channel.fromPath(params.cat_fasta))) | remap_segment
        
        remap_gtfs = remapped.map{ tup -> 
            tup[0] 
        }.collect() | cat_remap_gtfs
        remap_gtfs2 = remap_gtfs.combine(svgenes_out.map{ tup -> tup[2] }) | cat_remap_gtfs2
        remap_genes = remapped.map{ tup -> 
            tup[1]
        }.collect() | cat_remap_genes
        
        svgenes_gtf = svgenes_out.map{ tup ->
            def match = ( tup[2] =~ /([^\/]+).(gtf|gff|gff3)/)[0]
            return [match[0], match[1], tup[2] ]
        }
         
        svgenes2_out = gtf1.combine(svgenes_gtf).combine(Channel.fromPath(params.hg38_fasta)).combine(Channel.fromPath(params.cat_fasta)) | svgenes2
       
        genesrm.combine(svgenes_out.map{ tup -> [tup[0], tup[4]] }).combine(svgenes2_out.map{ tup -> 
            return [tup[2], tup[4]]
        }).combine(remap_genes) | finalize
    }
}

