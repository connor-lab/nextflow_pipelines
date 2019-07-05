#!/usr/bin/env nextflow

// INPUT READS
allFastq = "${params.dir}/Data/Intensities/BaseCalls/"

// INPUT RUN FOLDER
runDir = "${params.dir}"

// OUTPUT DIR
outDir = "${params.outdir}"

// MultiQC config file
multiQCConf = "${params.multiqcconf}"

// Centrifuge DB
centrifugeDB = "${params.centrifugedbpath}"

// SORT READS BY SIZE AND TAKE THE LARGEST ONE - MAKES SURE IT HAS AT LEAST ONE READ
allFilesBySize = file(allFastq).listFiles().sort{ it.size() }.reverse()
largestFile = allFilesBySize[0]

// THIS IS HORRIBLE BUT IT WORKS - CALL BASH SCRIPT TO GET RUNID
RunIDCmd = "getRunID ${largestFile}"
RunID = "${RunIDCmd}".execute().text.replaceAll(/\n/, "").toString()

// PROJECT ID LIST - these should map to shares and are defined in config file
projects = params.projectlist


// Setup input channels
// This takes all fastqs and filters so that only those defined in config file (projectlist) get processed
Channel.fromFilePairs( "${allFastq}/*_R{1,2}_001.fastq.gz" , flat: true)
       .filter { it[0].tokenize("-")[0].toLowerCase() in projects }
       //.into{ CopyRawFastq; MakeProjectFastq }
       .set{ InputReads }


// Channel for MiSeq run directory, needed for InterOp process
Channel.fromPath( runDir, type: 'dir', maxDepth: 1 )
       .into{ IlluminaRunDir ; backupRunDir }

// Channel for centrifuge database files, so that they are symlinked into WORKDIR
Channel
    .fromPath( centrifugeDB )
    .set{ CentrifugeRefs }

// Channel for NCBI taxIDs for each project
Channel
    .from( params.taxIDdict )
    .set { taxID }

// Channel for multiqc config file
Channel
    .fromPath( multiQCConf )
    .set{  multiQCConfYaml }

process IlluminaInteropStats {
    tag{ RunID }

    container "file:///${params.simgdir}/interop.simg"

    input:
    file(interopDir) from IlluminaRunDir

    output:
    file "*.summary.csv" into InterOpMultiQC

    script:
      """
      interop_summary --csv=1 ${interopDir} > ${RunID}.summary.csv
      """
}

backupHost = params.backupHost
backupUser = params.backupUser
backupPath = params.backupPath


process BackupRunDirectory {
    tag{RunID}

    input:
    file(runDir) from backupRunDir

    script:
    miseqID = "${runDir}".tokenize("_")[1]
    """
    rsync -avL ${runDir} ${backupUser}@${backupHost}:${backupPath}/${miseqID}/
    """
}

InputReads.map{ [ it[0] , it[0].tokenize("-")[0].toLowerCase() , it[1] , it[2] ]}.into{ inputReadsProject ; countProject }

countProject.map{ [ it[0] , it[1] ]}. countBy{ it[1] }.println()

process TrimReads {
    tag { dataset_id }

    container "file:///${params.simgdir}/trim_galore.simg"
   
    publishDir "${outDir}/${project}/${RunID}/qc/trimmed_reads", pattern: '*_val_{1,2}.fq.gz', mode: 'copy'
    publishDir "${outDir}/${project}/${RunID}/qc/fastqc", pattern: '*_fastqc.{zip,html}', mode: 'copy'
    publishDir "${outDir}/${project}/${RunID}/qc/trim_galore" , pattern: '*_trimming_report.txt', mode: 'copy'

    cpus 2

    input: 
    set dataset_id, project, file(forward), file(reverse) from inputReadsProject
 
    output:
    set dataset_id, project, file("*_val_1.fq.gz"), file("*_val_2.fq.gz") optional true into TrimmedReadsLength, TrimmedReadsQC, TrimmedReadsHIV, TrimmedReadsFLU, TrimmedReadsWCMTB, TrimmedReadsARG, TrimmedReadsDIGCD
    set project, file("*trimming_report.txt") optional true into TrimGaloreResults
    set project, file("*_fastqc.{zip,html}") optional true into TrimGaloreFastQCReports
    

    script:
      """
      if [[ \$(zcat ${forward} | head -n4 | wc -l) -eq 0 ]]; then
        exit 0
      else
        trim_galore --fastqc --paired $forward $reverse
      fi
      """
}

process MeanTrimmedReadLength {
    tag { dataset_id }

    container "file:///${params.simgdir}/seqtk.simg"

    input:
    set dataset_id, project, file(forward), file(reverse) from TrimmedReadsLength

    output:
    set project, dataset_id, stdout into ReadsLengthSummary

    script:
      """
      zcat $forward $reverse | seqtk fqchk - | head -n1 | cut -d ";" -f3 | cut -d " " -f3 | tr -d '\n'
      """
}


process Centrifuge {
    tag { dataset_id }

    container "file:///${params.simgdir}/centrifuge.simg"

    publishDir "${outDir}/${project}/${RunID}/qc/centrifuge", mode: 'copy'

    cpus 8

    //queue 'compute-hm'

    input:
    set dataset_id, project, file(forward), file(reverse) from TrimmedReadsQC
    file dbs from CentrifugeRefs.toList()

    output:
    set project, file("*_centrifugereport.tab") into CentrifugeReport
    set project, dataset_id, file("*_centrifugereport.tab") into CentrifugeSummary

    script:
      """
      centrifuge --mm -q -p ${task.cpus} -x phw -1 $forward -2 $reverse -S /dev/null --report-file ${dataset_id}_centrifugereport.tab
      """
}

process CentrifugeSummary {
    tag { dataset_id }

    container "file:///${params.simgdir}/taxonkit.simg"

    cpus 2

    input:
    set project, dataset_id, file(centrifuge_report), taxID from CentrifugeSummary.combine(taxID, by:0)

    output:
    set dataset_id, project, stdout into TotalBpCalc
    set project, taxID into MagnitudePrepare, TranslateTaxonomy

    script:
      """
      centrifuge-summary -a -i ${taxID} -r ${centrifuge_report} | sed 's/No matching reads found in this sample/0/g' | cut -f2 | tr -d '\n'
      """
}

process TranslateTaxonomy {
    tag { project }

    container "file:///${params.simgdir}/taxonkit.simg"

    cpus 2

    input:
    set project, taxID from TranslateTaxonomy.unique()

    output:
    set project, stdout into MagnitudeSummaryTaxName

    script:
      """
      taxonkit list --indent "" --show-name --ids ${taxID} | head -n1 | cut -d " " -f2- | tr -d '\n'
      """
}


process PrepareMagnitudeSummary {
    tag { project }

    executor 'local'

    input:
    set project, taxID from MagnitudePrepare.unique()

    output:
    set taxID, project into MagnitudeSummaryTaxID

    exec:
    def file = new File("${outDir}/${project}/${RunID}/qc/magnitude_summary.csv")
    if (file.exists()){
      file.delete();
    }
    file.append("Sample,Taxon,NCBI TaxID,Number of reads matching taxon,Mean read length,Total bp attributed to taxon\n")
}


process CalculateMagnitudeSummary {
    tag { dataset_id }

    executor 'local'

    input:
    set project, dataset_id, numberReads, taxID, taxName, avgReadLength from TotalBpCalc.combine(MagnitudeSummaryTaxID, by: 1).combine(MagnitudeSummaryTaxName, by: 0).combine(ReadsLengthSummary, by:[0,1])

    exec:
    matching = Float.parseFloat(numberReads) * Float.parseFloat(avgReadLength)
    matchingBP = matching.toInteger()
    def file = new File("${outDir}/${project}/${RunID}/qc/magnitude_summary.csv")
    file.append("${dataset_id},${taxName},${taxID},${numberReads},${avgReadLength},${matchingBP}\n")
}

process FlattenOutputs {
    tag { RunID }

    executor 'local'

    input:
    set project, trimoutputs from TrimGaloreResults.groupTuple()
    set project, fastqcoutputs from TrimGaloreFastQCReports.groupTuple()
    set project, centrifugeoutputs from CentrifugeReport.groupTuple()

    output:
    set project, trimflat, fastqcflat into MultiQCFlat
    set project, centrifugeflat into KronaFlat

    exec:
    trimflat = trimoutputs.flatten()
    fastqcflat = fastqcoutputs.flatten()
    centrifugeflat = centrifugeoutputs.flatten()
}

process MultiQC {
    tag { proctag }

    container "file:///${params.simgdir}/multiqc.simg"

    publishDir "${outDir}/${project}/${RunID}/qc", mode: 'copy'

    input:
    set project, file(trim), file(fastqc), file(interop), file(multiqcconf) from MultiQCFlat.combine(InterOpMultiQC).combine(multiQCConfYaml)

    output:
    file "multiqc_report.html"
    file "*_data"

    script:
    projectUpper = "${project}".toUpperCase()
    proctag = RunID + "-" + projectUpper
      """
      multiqc -m cutadapt -m fastqc -m interop -i "${projectUpper} ${RunID}" -n multiqc_report.html .
     """
}


process Krona {
    tag { proctag }

    container "file:///${params.simgdir}/kronatools.simg"

    publishDir "${outDir}/${project}/${RunID}/qc", mode: 'copy'

    input:
    set project, file(centrifuge) from KronaFlat

    output:
    file "centrifuge_report.html"

    script:
      projectUpper = "${project}".toUpperCase()
      proctag = RunID + "-" + projectUpper
      """
      ktImportTaxonomy -o centrifuge_report.html -m 5 -s 7 *_centrifugereport.tab
      """
}

// ** ## HIV PIPELINE ## ** //

// Setup HIV trimmed reads channel
HIVTrimmedReads = TrimmedReadsHIV.filter { it[1] == 'hiv' }

// Get today's year and month
Date date = new Date()
YearMonth = date.format( "yyyy/MM-MMM" )


// Make reports directory if it doesn't exist
ReportsDir = file("/mnt/datastore/hiv/reports/${YearMonth}")
mkdirResult = ReportsDir.mkdirs()
println mkdirResult ? "Made new reports directory" : "Cannot create directory: $ReportsDir"

// Setup references
HIVComp = Channel.fromPath( "${params.subref}" )
HIVHXB2 = Channel.fromPath( "${params.HXB2ref}" )

// Setup init_dir for shiver
ShiverInit = Channel.fromPath( "${params.shiverinit}", type: 'dir', maxDepth: 1)
ShiverConf = Channel.fromPath( "${params.shiverconf}" )

// Minor variant frequency list
MinVarFreq = Channel.from(params.minvarfreq)

process CleanHIVReads {
    tag { dataset_id }
  
    container "file:///${params.simgdir}/minimap2.simg"

    cpus 4

    input:
    set dataset_id, project, file(forward), file(reverse), file(ref) from HIVTrimmedReads.combine(HIVComp)

    output:
    set dataset_id, file("${dataset_id}.clean_1.fq.gz"), file("${dataset_id}.clean_2.fq.gz") into HIVCleanReadsAssembly, HIVCleanReadsPolishing, HIVCleanReadsVariantCalling

    script:
    """
    minimap2 -t ${task.cpus} -x sr -a $ref $forward $reverse | samtools view -F 4 -@ 2 -b > clean.bam
    picard SamToFastq VALIDATION_STRINGENCY=LENIENT I=clean.bam F=${dataset_id}.clean_1.fq.gz F2=${dataset_id}.clean_2.fq.gz
    """
}

process SampleHIVReads {
    tag { dataset_id }

    container "file:///${params.simgdir}/bbtools.simg"

    publishDir "${params.hivdir}/${RunID}/analysis/01-clean_subsampled_reads", pattern: '*.clean.sampled_{1,2}.fq.gz', mode: 'copy'

    input:
    set dataset_id, file(forward), file(reverse) from HIVCleanReadsAssembly

    output:
    set dataset_id, file("${dataset_id}.clean.sampled_1.fq.gz"), file("${dataset_id}.clean.sampled_2.fq.gz") into HIVSampledReadsAssembly

    script:
    """
    reformat.sh samplebasestarget=${params.samplebases} in=$forward in2=$reverse out=formatted_1.fastq.gz out2=formatted_2.fastq.gz
    shuffle.sh -Xmx2g in=formatted_1.fastq.gz in2=formatted_2.fastq.gz out=${dataset_id}.clean.sampled_1.fq.gz out2=${dataset_id}.clean.sampled_2.fq.gz
    """
}


process AssembleHIVReads {
    tag { dataset_id }

    publishDir "/mnt/datastore/hiv/${RunID}/analysis/02-assembly", pattern: "${dataset_id}.iva.fa", mode: 'copy'
 
    container "file:///${params.simgdir}/iva.simg"

    cpus 4

    input:
    set dataset_id, file(forward), file(reverse) from HIVSampledReadsAssembly

    output:
    set dataset_id, file("${dataset_id}.iva.fa") optional true into HIVIVAAssembly
    set dataset_id, file("*.assembly.failed") optional true into HIVIVAFail

    script:
    """
    if iva --threads ${task.cpus} -f $forward -r $reverse iva_assembly &> ${dataset_id}.assembly.failed ; then
      mv iva_assembly/contigs.fasta ${dataset_id}.iva.fa
    else
      sed -i "1s/^/${dataset_id}\t/" ${dataset_id}.assembly.failed
    fi
    """
}


process CollectFailedHIVAssemblies {
    tag { RunID }

    publishDir "${params.hivdir}/${RunID}/analysis/", pattern: "IVA.assembly.failed.txt", mode: 'copy'

    input:
    file "*" from HIVIVAFail.collect()

    output:
    file "IVA.assembly.failed.txt"

    script:
    """
    cat *.failed > IVA.assembly.failed.txt
    """
}


if(params.shiver == 'false'){
process OrderHIVContigs {
    tag { dataset_id }

    container "file:///${params.simgdir}/assembly_improvement.simg"

    input:
    set dataset_id, file(assembly), file(ref) from HIVIVAAssembly.combine(HIVHXB2)

    output:
    set dataset_id, file("${dataset_id}.ordered.fa") into HIVIVAAssemblyOrdered

    script:
    """
    order_contigs_with_abacas -a $assembly -c $ref
    mv ${dataset_id}.iva.fa.scaffolded.filtered ${dataset_id}.ordered.fa
    sed -i "s/>/>${dataset_id}_ordered_/g" ${dataset_id}.ordered.fa
    """
}

process GapfillHIVContigs {
    tag { dataset_id }

    publishDir "${params.hivdir}/${RunID}/analysis/02-assembly", pattern: "${dataset_id}.polished.fa", mode: 'copy'

    container "file:///${params.simgdir}/assembly_improvement.simg"

    cpus 4

    input:
    set dataset_id, file(assembly), file(forward), file(reverse) from HIVIVAAssemblyOrdered.join(HIVCleanReadsPolishing)
    
    output:
    set dataset_id, file("${dataset_id}.polished.fa") into HIVAssemblyBAM, HIVAssemblyVariants

    script:
    // NEED TO TOUCH FILES OR GAP2SEQ FALLS OVER ON CONTIGS WITH NO GAP
    """
    touch tmp.gaps
    touch tmp.fill
    Gap2Seq.sh -nb-cores ${task.cpus} -k 51 -scaffolds $assembly -filled ${dataset_id}.polished.fa -reads ${forward},${reverse}
    sed -i "s/ordered/polished/g" ${dataset_id}.polished.fa
    """
}}
else {

process HIVShiver {
    tag { dataset_id }

    container "file:///${params.simgdir}/shiver.simg"

    publishDir "${params.hivdir}/${RunID}/analysis/02-assembly", pattern: "${dataset_id}.shiver.fa", mode: 'copy'
   
    input:
    set dataset_id, file(assembly), file(forward), file(reverse), file(shiverconf), file(shiverinit) from HIVIVAAssembly.join(HIVCleanReadsPolishing).combine(ShiverConf).combine(ShiverInit)

    output:
    set dataset_id, file("${dataset_id}.shiver.fa") into HIVAssemblyBAM, HIVAssemblyVariants

    script:
    """
    shiver_align_contigs.sh ${shiverinit} ${shiverconf} ${assembly} ${dataset_id}
    if [ -f ${dataset_id}_cut_wRefs.fasta ]; then
      shiver_map_reads.sh ${shiverinit} ${shiverconf} ${assembly} ${dataset_id} ${dataset_id}.blast ${dataset_id}_cut_wRefs.fasta ${forward} ${reverse}
    else
      shiver_map_reads.sh ${shiverinit} ${shiverconf} ${assembly} ${dataset_id} ${dataset_id}.blast ${dataset_id}_raw_wRefs.fasta ${forward} ${reverse}
    fi
    seqtk seq -l0 ${dataset_id}_remap_consensus_MinCov_15_30.fasta | head -n2 | sed '/>/!s/-//g' | sed 's/\\?/N/g' | sed 's/_remap_consensus//g' | seqtk seq -l80 > ${dataset_id}.shiver.fa
    """
}
}


process HIVMappingVariantCalling {
    tag { dataset_id }

    cpus 4

    container "file:///${params.simgdir}/minimap2.simg"

    input:
    set dataset_id, file(forward), file(reverse), file(assembly) from HIVCleanReadsVariantCalling.join(HIVAssemblyBAM)

    output:
    set dataset_id, file("${dataset_id}.nodups.bam") into HIVMappingNoDupsBAM

    script:
    """
    samtools faidx $assembly
    minimap2 -t ${task.cpus} -x sr -a $assembly $forward $reverse | samtools view -@ 2 -b | samtools sort -@ 2 -o ${dataset_id}.variants.bam
    picard MarkDuplicates METRICS_FILE=duplicates_metrics.txt REMOVE_DUPLICATES=true I=${dataset_id}.variants.bam O=${dataset_id}.nodups.bam
    """
}


if(params.variantstrategy.toLowerCase() == 'varscan'){
process HIVVariantCallingVarScan {
    tag { dataset_id }

    container "file:///${params.simgdir}/variant_calling.simg"

    publishDir "${params.hivdir}/${RunID}/analysis/03-call_variants/fasta/minor_variants", pattern: "${dataset_id}.${minvarfreq}.minor.fa", mode: 'copy'
    publishDir "${params.hivdir}/${RunID}/analysis/03-call_variants/fasta/IUPAC", pattern: "${dataset_id}.${minvarfreq}.iupac.consensus.fa", mode: 'copy'
    publishDir "${params.hivdir}/${RunID}/analysis/03-call_variants/vcf", pattern: "${dataset_id}.${minvarfreq}.consensus.vcf", mode: 'copy'

    input:
    set dataset_id, file(nodupsbam), file(assembly), minvarfreq from HIVMappingNoDupsBAM.join(HIVAssemblyVariants).combine(MinVarFreq)

    output:
    set dataset_id, minvarfreq, file("*.${minvarfreq}.iupac.consensus.fa") into HIVAssemblyWithVariants
    file "${dataset_id}.${minvarfreq}.minor.fa"
    file "${dataset_id}.${minvarfreq}.consensus.vcf"

    script:
    """
    samtools mpileup --max-depth 100000 --redo-BAQ --min-MQ 17 --min-BQ 20 --output ${dataset_id}.mpileup --fasta-ref ${assembly} ${nodupsbam}
    java -Xmx17G -jar /usr/local/bin/varscan.jar mpileup2cns ${dataset_id}.mpileup --min-var-freq ${minvarfreq} --p-value 95e-02 --min-coverage 100 --output-vcf 1 > ${dataset_id}.varscan.cns.vcf
    bgzip ${dataset_id}.varscan.cns.vcf
    tabix -p vcf ${dataset_id}.varscan.cns.vcf.gz
    bcftools view -i'FILTER="PASS"' -Oz -o ${dataset_id}.varscan.cns.filtered.vcf.gz ${dataset_id}.varscan.cns.vcf.gz
    zcat ${dataset_id}.varscan.cns.filtered.vcf.gz > ${dataset_id}.${minvarfreq}.consensus.vcf
    tabix -p vcf ${dataset_id}.varscan.cns.filtered.vcf.gz
    bcftools consensus -f $assembly ${dataset_id}.varscan.cns.filtered.vcf.gz --output ${dataset_id}.${minvarfreq}.minor.fa
    bcftools consensus -I -f $assembly ${dataset_id}.varscan.cns.filtered.vcf.gz --output ${dataset_id}.${minvarfreq}.iupac.consensus.fa
    sed -i 's/polished/consensus-minor/g' ${dataset_id}.${minvarfreq}.minor.fa
    sed -i 's/polished/consensus-iupac/g' ${dataset_id}.${minvarfreq}.iupac.consensus.fa
    sed -i '/^>/ s/\$/ [Variant caller: ${params.variantstrategy}] [Minor variant bases IUPAC] [Variant frequency: ${minvarfreq}] /' ${dataset_id}.${minvarfreq}.iupac.consensus.fa
    sed -i '/^>/ s/\$/ [Variant caller: ${params.variantstrategy}] [Minor variants bases ONLY] [Variant frequency: ${minvarfreq}] /' ${dataset_id}.${minvarfreq}.minor.fa
    """
}
}


if(params.variantstrategy == 'bcftools'){
process HIVVariantCallingBCFtools {
    tag { dataset_id }

    container "file:///${params.simgdir}/variant_calling.simg"

    publishDir "${params.hivdir}/${RunID}/analysis/03-call_variants/fasta/minor_variants", pattern: "${dataset_id}.${minvarfreq}.minor.fa", mode: 'copy'
    publishDir "${params.hivdir}/${RunID}/analysis/03-call_variants/fasta/IUPAC", pattern: "${dataset_id}.${minvarfreq}.iupac.consensus.fa", mode: 'copy'
    publishDir "${params.hivdir}/${RunID}/analysis/03-call_variants/vcf", pattern: "${dataset_id}.${minvarfreq}.consensus.vcf", mode: 'copy'


    input:
    set dataset_id, file(nodupsbam), file(assembly), minvarfreq from HIVMappingNoDupsBAM.join(HIVAssemblyVariants).combine(MinVarFreq)

    output:
    set dataset_id, minvarfreq, file("*.${minvarfreq}.iupac.consensus.fa") into HIVAssemblyWithVariants
    file "${dataset_id}.${minvarfreq}.minor.fa"
    file "${dataset_id}.${minvarfreq}.consensus.vcf"

    script:
    """
    samtools faidx $assembly
    bcftools mpileup --redo-BAQ --min-MQ 20 -Ou -f $assembly ${dataset_id}.nodups.bam | bcftools call --ploidy 1 -mv -Ov | bcftools view -q ${minvarfreq}:nref -Oz -o ${dataset_id}.variants.vcf.gz
    zcat ${dataset_is}.variants.vcf.gz > ${dataset_id}.${minvarfreq}.consensus.vcf
    tabix ${dataset_id}.variants.vcf.gz
    bcftools consensus -f $assembly ${dataset_id}.variants.vcf.gz --output ${dataset_id}.${minvarfreq}.minor.fa
    bcftools consensus -I -f $assembly ${dataset_id}.variants.vcf.gz --output ${dataset_id}.${minvarfreq}.iupac.consensus.fa
    sed -i 's/polished/consensus-minor/g' ${dataset_id}.${minvarfreq}.minor.fa
    sed -i 's/polished/consensus-iupac/g' ${dataset_id}.${minvarfreq}.iupac.consensus.fa
    sed -i '/^>/ s/\$/ [Variant caller: ${params.variantstrategy}] [Minor variant bases IUPAC] [Variant frequency: ${minvarfreq}] /' ${dataset_id}.${minvarfreq}.iupac.consensus.fa
    sed -i '/^>/ s/\$/ [Variant caller: ${params.variantstrategy}] [Minor variants bases ONLY] [Variant frequency: ${minvarfreq}] /' ${dataset_id}.${minvarfreq}.minor.fa
    """
}
}

if(params.variantstrategy == 'lofreq'){
process HIVVariantCallingLoFreq {
    tag { dataset_id }

    container "file:///${params.simgdir}/variant_calling.simg"

    publishDir "${params.hivdir}/${RunID}/analysis/03-call_variants/fasta/minor_variants", pattern: "${dataset_id}.${minvarfreq}.minor.fa", mode: 'copy'
    publishDir "${params.hivdir}/${RunID}/analysis/03-call_variants/fasta/IUPAC", pattern: "${dataset_id}.${minvarfreq}.iupac.consensus.fa", mode: 'copy'
    publishDir "${params.hivdir}/${RunID}/analysis/03-call_variants/vcf", pattern: "${dataset_id}.${minvarfreq}.consensus.vcf", mode: 'copy'

    input:
    set dataset_id, file(nodupsbam), file(assembly), minvarfreq from HIVMappingNoDupsBAM.join(HIVAssemblyVariants).combine(MinVarFreq)

    output:
    set dataset_id, minvarfreq, file("*.${minvarfreq}.iupac.consensus.fa") into HIVAssemblyWithVariants
    file "${dataset_id}.${minvarfreq}.minor.fa"
    file "${dataset_id}.${minvarfreq}.consensus.vcf"

    script:
    """
    samtools faidx $assembly
    lofreq call -C 100 --call-indels -D -f $assembly -o - $nodupsbam | bcftools view -i'AF>${minvarfreq}' -Oz -o ${dataset_id}.variants.vcf.gz
    zcat ${dataset_id}.variants.vcf.gz > ${dataset_id}.${minvarfreq}.consensus.vcf
    tabix ${dataset_id}.variants.vcf.gz
    bcftools consensus -f $assembly ${dataset_id}.variants.vcf.gz --output ${dataset_id}.${minvarfreq}.minor.fa
    bcftools consensus -I -f $assembly ${dataset_id}.variants.vcf.gz --output ${dataset_id}.${minvarfreq}.iupac.consensus.fa
    sed -i 's/polished/consensus-minor/g' ${dataset_id}.${minvarfreq}.minor.fa
    sed -i 's/polished/consensus-iupac/g' ${dataset_id}.${minvarfreq}.iupac.consensus.fa
    sed -i '/^>/ s/\$/ [Variant caller: ${params.variantstrategy}] [Minor variant bases IUPAC] [Variant frequency: ${minvarfreq}] /' ${dataset_id}.${minvarfreq}.iupac.consensus.fa
    sed -i '/^>/ s/\$/ [Variant caller: ${params.variantstrategy}] [Minor variants bases ONLY] [Variant frequency: ${minvarfreq}] /' ${dataset_id}.${minvarfreq}.minor.fa
    """
}
}

HIVAssemblyWithSelectedMinVarFreq = HIVAssemblyWithVariants.filter{ it[1] == params.selectedminvarfreq }

process HIVMakeResistanceReport {
    tag { dataset_id }

    container "file:///${params.simgdir}/sierrapy.simg"

    publishDir "${params.hivdir}/${RunID}/analysis/04-call_resistance", pattern: "${dataset_id}.json", mode: 'copy'
    publishDir "${params.hivdir}/${RunID}/analysis/05-generate_report", pattern: "${dataset_id}.rtf", mode: 'copy'
    publishDir "${params.hivdir}/reports/${YearMonth}", pattern: "${dataset_id}.rtf", mode: 'copy'

    //queue 'internet'

    input:
    set dataset_id, minvarfreq, file(variantassembly) from HIVAssemblyWithSelectedMinVarFreq
//    set dataset_id, file(stanfordjson) from HIVDBReportJson

    output:
    file("${dataset_id}.json")
    file("${dataset_id}.rtf")

    script:
    episodenumber = dataset_id.split('-').last()
    """
    buildreport.pl -i ${variantassembly} -n ${dataset_id} -l PHW_Cardiff
    """
}

// ** ## FLU PIPELINE ## ** //

// Setup FLU trimmed reads channel
FLUTrimmedReads = TrimmedReadsFLU.filter { it[1] == 'flu' }

// Setup init_dir for shiver
FLUShiverInit = Channel.fromPath( "${params.flushiverinitroot}", type: 'dir')
FLUShiverConf = Channel.fromPath( "${params.flushiverconf}" )

// Setup references
FLURef = Channel.from(['NA', file(params.flurefNA)],
                      ['HA', file(params.flurefHA)],
                      ['M1', file(params.flurefM1)],
                      ['PB1', file(params.flurefPB1)],
                      ['PB2', file(params.flurefPB2)],
                      ['NP', file(params.flurefNP)],
                      ['NS1', file(params.flurefNS1)],
                      ['PA', file(params.flurefPA)])


process SeparateFLUSegmentReads {
    tag { proctag }

    container "file:///${params.simgdir}/minimap2.simg"

    cpus 4

    input:
    set dataset_id, project, file(forward), file(reverse), segment, file(ref) from FLUTrimmedReads.combine(FLURef)

    output:
    set dataset_id, segment, file("${dataset_id}.${segment}_1.fq.gz"), file("${dataset_id}.${segment}_2.fq.gz") into FLUCleanReadsAssembly,FLUCleanReadsShiver

    script:
    proctag = dataset_id + "-" + segment
    """
    minimap2 -t ${task.cpus} -x sr -a $ref $forward $reverse | samtools view -F 4 -@ 2 -b > mapped.bam
    picard SamToFastq VALIDATION_STRINGENCY=LENIENT I=mapped.bam F=${dataset_id}.${segment}_1.fq.gz F2=${dataset_id}.${segment}_2.fq.gz
    """
}


process AssembleFLUReads {
    tag { proctag }

    container "file:///${params.simgdir}/iva.simg"

    cpus 4

    input:
    set dataset_id, segment, file(forward), file(reverse) from FLUCleanReadsAssembly

    output:
    set dataset_id, segment, file("${dataset_id}.${segment}.iva.fa") optional true into FLUIVAAssembly
    set dataset_id, file("*.assembly.failed") optional true into FLUIVAFail

    script:
    proctag = dataset_id + "-" + segment
    """
    if iva --threads ${task.cpus} -f $forward -r $reverse iva_assembly &> ${dataset_id}.${segment}.assembly.failed ; then
      mv iva_assembly/contigs.fasta ${dataset_id}.${segment}.iva.fa
      sed -i "s/>/>${dataset_id}\\.${segment}\\./g" ${dataset_id}.${segment}.iva.fa
    else
      sed -i "1s/^/${dataset_id}\\.${segment}\t/" ${dataset_id}.${segment}.assembly.failed
    fi
    """
}


fluAssemblyShiver = Channel.create()

fluAssemblyNotShiver = Channel.create()

FLUIVAAssembly.choice( fluAssemblyShiver, fluAssemblyNotShiver){ it[1] in params.shiversegs ? 0 : 1 }


process shiverFLU {
    tag { proctag }

    container "file:///${params.simgdir}/shiver.simg"

    cpus 4

    input:
    set dataset_id, segment, file(assembly), file(forward), file(reverse), file(shiverconf), file(shiverinit) from fluAssemblyShiver.combine(FLUCleanReadsShiver, by: [ 0 , 1 ]).combine(FLUShiverConf).combine(FLUShiverInit)

    output:
    set dataset_id, file("${dataset_id}.${segment}.shiver.fa") into fluAssemblyReconstitute

    script:
    proctag = dataset_id + "-" + segment
    """
    shiver_align_contigs.sh ${shiverinit}/${segment} ${shiverconf} ${assembly} ${dataset_id}
    if [ -f ${dataset_id}_cut_wRefs.fasta ]; then
      shiver_map_reads.sh ${shiverinit}/${segment} ${shiverconf} ${assembly} ${dataset_id} ${dataset_id}.blast ${dataset_id}_cut_wRefs.fasta ${forward} ${reverse}
    else
      shiver_map_reads.sh ${shiverinit}/${segment} ${shiverconf} ${assembly} ${dataset_id} ${dataset_id}.blast ${dataset_id}_raw_wRefs.fasta ${forward} ${reverse}
    fi
    seqtk seq -l0 ${dataset_id}_remap_consensus_MinCov_15_30.fasta | head -n2 | sed '/>/!s/-//g' | sed 's/\\?/N/g' | sed "s/_remap_consensus/\\.${segment}/g" > ${dataset_id}.${segment}.shiver.fa
    """
}


process removeSegmentFLU {
    tag{ proctag }

    executor 'local'

    stageInMode 'copy'

    input:
    set dataset_id, segment, file(assembly) from fluAssemblyNotShiver

    output:
    set dataset_id, file("${dataset_id}.${segment}.iva.fa") into fluAssemblyNotShiverClean

    script:
    proctag = dataset_id + "-" + segment
    """
    """
}


process ReconstituteFLUGenome {
    tag { dataset_id }

    container "file:///${params.simgdir}/seqtk.simg"

    publishDir "${params.fludir}/${RunID}/analysis/assembly", pattern: "${dataset_id}.fasta", mode: 'copy'

    input:
    set dataset_id, file('*') from fluAssemblyReconstitute.concat(fluAssemblyNotShiverClean).groupTuple()

    output:
    file "${dataset_id}.fasta"

    script:
    """
    cat *.fa | seqtk seq -l0 > ${dataset_id}.fasta
    """
}


process CollectFailedFLUAssemblies {
    tag { RunID }

    publishDir "${params.fludir}/${RunID}/analysis/", pattern: "IVA.assembly.failed.txt", mode: 'copy'

    input:
    file "*" from FLUIVAFail.collect()

    output:
    file "IVA.assembly.failed.txt"

    script:
    """
    cat *.failed > IVA.assembly.failed.txt
    """
}

// ** ## WCM PIPELINE ## ** //

// Setup WCM trimmed reads channel
WCMTBTrimmedReads = TrimmedReadsWCMTB.filter { it[1] == 'wcmtb' || it[1] == 'wcmid' }

WCMTBTrimmedReads.into{ WCMTBTrimmedReadsMD5; WCMTrimmedReadsTyping; WCMMakeUploadDir; WCMTrimmedReadsKraken; WCMTrimmedReadsShovill }

mykrobepanel = params.mykrobepanel

WCMKrakenDB = params.wcmkrakendbdir

Channel
    .fromPath( WCMKrakenDB )
    .set{ WCMKrakenRefs }


Date day = new Date()
YearMonthDay = day.format( "yyyy-MM-dd" )

shortRunID = RunID.tokenize("_")[0, 1].join("_")

sshkey = "${params.wcmtbkey}"
endserver = "${params.wcmtbserver}"
user = "${params.wcmtbuser}"
endpath = "${params.wcmtbpath}/${YearMonthDay}-${shortRunID}"


process assembleShovillWCM {
    tag { dataset_id }

    errorStrategy 'ignore'

    container "file:///${params.simgdir}/shovill.simg"

    cpus 5

    memory '8 GB'

    input:
    set dataset_id, project, file(forward), file(reverse) from WCMTrimmedReadsShovill

    output:
    set dataset_id, project, file("${dataset_id}.fasta") optional true into WCMProkkaAssembly
    set project, file("${dataset_id}.fasta") optional true into WCMquast
    set dataset_id, project, file("${dataset_id}.fasta") into WCMkrakenAssembly

    script:
    """
    shovill --cpus ${task.cpus} --R1 ${forward} --R2 ${reverse} --minlen 500 --outdir shovill
    mv shovill/contigs.fa ${dataset_id}.fasta
    """
}

process quastWCM {
    tag { project }

    publishDir "${outDir}/${project}/${RunID}/analysis", pattern: "*_assembly_summary.csv", mode: 'copy'

    container "file:///${params.simgdir}/quast.simg"

    cpus 4

    input:
    set project, file(assembly) from WCMquast.groupTuple(by: 0)

    output:
    file("*_assembly_summary.csv")

    script:
    """
    quast.py -t ${task.cpus} -o . --contig-thresholds 0 --no-html --no-plots *.fasta
    sed 's/\t/,/g' transposed_report.tsv > ${shortRunID}_assembly_summary.csv
    """
}


process annotateProkkaWCM {
    tag { dataset_id }

    publishDir "${outDir}/${project}/${RunID}/analysis/annotation/", pattern: "${dataset_id}/${dataset_id}.*", mode: 'copy'

    container "file:///${params.simgdir}/prokka.simg"

    cpus 5

    input:
    set dataset_id, project, file(assembly) from WCMProkkaAssembly.filter{ it[2].size()>1000 }

    output:
    file("${dataset_id}/${dataset_id}.*")
    set project, file("${dataset_id}.fna") into WCMresAssembly, WCMvirAssembly

    script:
    locusTag = dataset_id.tokenize("_")[0].tokenize("-")[1]
    prefix = project.toUpperCase()
    """
    prokka --outdir ${dataset_id} --locustag ${prefix}_${locusTag} --prefix ${dataset_id} --centre PHW --cpus ${task.cpus} --compliant ${assembly}
    cp ${dataset_id}/${dataset_id}.fna ${dataset_id}.fna
    """
}

process callResistanceWCM {
    tag { project }

    publishDir "${outDir}/${project}/${RunID}/analysis", pattern: "*resistance_genes.csv", mode: 'copy'

    container "file:///${params.simgdir}/abricate.simg"

    input:
    set project, file(assembly) from WCMresAssembly.groupTuple(by: 0)

    output:
    file "*_resistance_genes.csv"

    script:
    """
    abricate --db ${params.resistancedb} --csv *.fna > ${shortRunID}_resistance_genes.csv
    """
}

process callVirulenceWCM {
    tag { project }

    publishDir "${outDir}/${project}/${RunID}/analysis", pattern: "*_virulence_genes.csv", mode: 'copy'

    container "file:///${params.simgdir}/abricate.simg"

    input:
    set project, file(assembly) from WCMvirAssembly.groupTuple(by: 0)

    output:
    file "*_virulence_genes.csv"

    script:
    """
    abricate --db ${params.virulencedb} --csv *.fna > ${shortRunID}_virulence_genes.csv
    """
}

process krakenWCM {
    tag { dataset_id }

    publishDir "${outDir}/${project}/${RunID}/analysis/kraken", pattern: "${dataset_id}_krakenreport.txt", mode: 'copy'

    container "file:///${params.simgdir}/kraken2.simg"

    cpus 5

    input:
    set dataset_id, project, file(assembly) from WCMkrakenAssembly
    file dbs from WCMKrakenRefs.toList()

    output:
    file("*_krakenreport.txt")
    set project, file("${dataset_id}.tab") into WCMKrakenReport

    script:
    """
    sed -i 's/ /|/g' ${assembly}
    kraken2 --db . --threads ${task.cpus} --output ${dataset_id}_kraken.txt --report ${dataset_id}_krakenreport.txt ${assembly}
    sed 's/|/\\t/g' ${dataset_id}_kraken.txt | sed 's/len=//g' | cut -f2,3,7 > ${dataset_id}.tab
    #cut -f2,3 ${dataset_id}_kraken.tab > ${dataset_id}.tab
    """
}

process kronaWCM {
    tag { proctag }

    container "file:///${params.simgdir}/kronatools.simg"

    publishDir "${outDir}/${project}/${RunID}/analysis/", pattern: "${shortRunID}_wcm-kraken_report.html", mode: 'copy'

    input:
    set project, file(centrifuge) from WCMKrakenReport.groupTuple(by: 0)

    output:
    file "${shortRunID}_wcm-kraken_report.html"

    script:
      projectUpper = "${project}".toUpperCase()
      proctag = RunID + "-" + projectUpper
      """
      ktImportTaxonomy -m 2 -t 3 -q 1 -o ${shortRunID}_wcm-kraken_report.html *.tab
      """
}


process typingMykrobeWCM {
    tag { dataset_id }

    publishDir "${outDir}/${project}/${RunID}/analysis/typing/json", pattern: "${dataset_id}.json", mode: 'copy'
    publishDir "${outDir}/${project}/${RunID}/analysis/typing/csv", pattern: "${dataset_id}.csv", mode: 'copy'

    container "file:///${params.simgdir}/mykrobe-atlas.simg"

    errorStrategy 'ignore'

    cpus 1

    input:
    set dataset_id, project, file(forward), file(reverse) from WCMTrimmedReadsTyping

    output:
    set dataset_id, project, file("${dataset_id}.json") into WCMTypingResistanceReport
    file("${dataset_id}.csv")

    script:
    """
    mykrobe predict ${dataset_id} tb --panel ${mykrobepanel} --seq *.fq.gz --format json --output ${dataset_id}.json
    json_to_tsv.py ${dataset_id}.json | sed 's/\\t/,/g' > ${dataset_id}.csv
    """
}


if(params.wcmtransferupload == 'true'){
process makeUploadDir {
    tag { RunID }

    executor 'local'

    input:
    set dataset_id, project, file(forward), file(reverse) from WCMMakeUploadDir.first()

    exec:
    createfqdirectory ="ssh -i ${sshkey} ${user}@${endserver} mkdir -p ${endpath}/fastq"
    createfqdirectory.execute()

    createrepdirectory ="ssh -i ${sshkey} ${user}@${endserver} mkdir -p ${endpath}/reports"
    createrepdirectory.execute()
}


process WCMGenerateMD5 {
    tag { dataset_id }

    cpus 1

    input:
    set dataset_id, project, file(forward), file(reverse) from WCMTBTrimmedReadsMD5.filter{ it[0].toUpperCase() =~ /NTC/ ? false : true }

    output:
    set file("*_1.qc.md5.txt"), file("*_2.qc.md5.txt") into MD5s
    set dataset_id, file("*_1.fq.gz"), file("*_2.fq.gz") into UploadReads

    script:
       accession = dataset_id.tokenize("_")[0].tokenize("-")[1]
       episode = dataset_id.tokenize("_")[0].tokenize("-")[2]
       shortRunID = RunID.tokenize("_")[0, 1].join("_")
       if (accession == "XXX") {
           prefix = episode + "!" + episode
       } else if (accession == "CON") {
           prefix = episode + "!" + episode
       } else {
           prefix = accession + "!" + episode
       }
       """
       mv $forward ${prefix}-${shortRunID}_1.fq.gz
       mv $reverse ${prefix}-${shortRunID}_2.fq.gz
       md5sum ${prefix}-${shortRunID}_1.fq.gz > ${prefix}-${shortRunID}_1.qc.md5.txt
       md5sum ${prefix}-${shortRunID}_2.fq.gz > ${prefix}-${shortRunID}_2.qc.md5.txt
       """
}

process WCMUploadFiles {
    tag { dataset_id }

    cpus 1

    //queue 'internet'

    maxForks 3
    errorStrategy 'retry'
    maxRetries 5
    maxErrors 25

    input:
    set dataset_id, file(forward), file(reverse) from UploadReads
    set file(forwardmd5), file(reversemd5) from MD5s

    script:
      """
      scp -i ${sshkey} $forward ${user}@${endserver}:${endpath}/fastq
      scp -i ${sshkey} $forwardmd5 ${user}@${endserver}:${endpath}/fastq
      scp -i ${sshkey} $reverse ${user}@${endserver}:${endpath}/fastq
      scp -i ${sshkey} $reversemd5 ${user}@${endserver}:${endpath}/fastq
      ssh -i ${sshkey} ${user}@${endserver} \"cd ${endpath}/fastq; md5sum -c ${forwardmd5}\"
      ssh -i ${sshkey} ${user}@${endserver} \"cd ${endpath}/fastq; md5sum -c ${reversemd5}\"
      """
}


}

// ** ## ARGENT PIPELINE ## ** //


// Setup ARG trimmed reads channel
ARGTrimmedReads = TrimmedReadsARG.filter { it[1] ==~ 'argid' || it[1] ==~ 'argab' }


process assembleShovillARGENT {
    tag { dataset_id }

    errorStrategy 'ignore'

//  publishDir "${outDir}/${project}/${RunID}/analysis/assembly", pattern: "${dataset_id}.fasta", mode: 'copy'

    container "file:///${params.simgdir}/shovill.simg"

    cpus 5

    memory '8 GB'

    input:
    set dataset_id, project, file(forward), file(reverse) from ARGTrimmedReads

    output:
    set dataset_id, project, file("${dataset_id}.fasta") optional true into ProkkaAssembly
    set project, file("${dataset_id}.fasta") optional true into quastARGENT

    script:
    """
    shovill --cpus ${task.cpus} --R1 ${forward} --R2 ${reverse} --minlen 500 --outdir shovill
    mv shovill/contigs.fa ${dataset_id}.fasta
    """
}

process quastARGENT {
    tag { project }

    publishDir "${outDir}/${project}/${RunID}/analysis", pattern: "*_assembly_summary.csv", mode: 'copy'

    container "file:///${params.simgdir}/quast.simg"

    cpus 4

    input:
    set project, file(assembly) from quastARGENT.groupTuple(by: 0)

    output:
    file("*_assembly_summary.csv")

    script:
    """
    quast.py -t ${task.cpus} -o . --contig-thresholds 0 --no-html --no-plots *.fasta
    sed 's/\t/,/g' transposed_report.tsv > ${RunID}_assembly_summary.csv
    """
}


process annotateProkkaARGENT {
    tag { dataset_id }

    publishDir "${outDir}/${project}/${RunID}/analysis/annotation/", pattern: "${dataset_id}/${dataset_id}.*", mode: 'copy'

    container "file:///${params.simgdir}/prokka.simg"

    cpus 5

    input:
    set dataset_id, project, file(assembly) from ProkkaAssembly.filter{ it[2].size()>1000 }

    output:
    file("${dataset_id}/${dataset_id}.*")
    set project, file("${dataset_id}.fna") into MLSTAssemblyARGENT, resAssemblyARGENT, virAssemblyARGENT

    script:
    locusTag = dataset_id.tokenize("_")[0].tokenize("-")[1]
    prefix = project.toUpperCase()
    """
    prokka --outdir ${dataset_id} --locustag ${prefix}_${locusTag} --prefix ${dataset_id} --centre PHW --cpus ${task.cpus} --compliant ${assembly}
    cp ${dataset_id}/${dataset_id}.fna ${dataset_id}.fna
    """
}

process callMLSTARGENT {
    tag { project }

    publishDir "${outDir}/${project}/${RunID}/analysis/", pattern: "*_MLST.csv", mode: 'copy'

    container "file:///${params.simgdir}/mlst.simg"

    input:
    set project, file(assembly) from MLSTAssemblyARGENT.groupTuple(by: 0)

    output:
    file "*_MLST.csv"

    script:
    """
    mlst --nopath --csv *.fna > ${RunID}_MLST.csv
    """
}

process callResistanceARGENT {
    tag { project }

    publishDir "${outDir}/${project}/${RunID}/analysis", pattern: "*resistance_genes.csv", mode: 'copy'

    container "file:///${params.simgdir}/abricate.simg"

    input:
    set project, file(assembly) from resAssemblyARGENT.groupTuple(by: 0)

    output:
    file "*_resistance_genes.csv"

    script:
    """
    abricate --db ${params.resistancedb} --csv *.fna > ${RunID}_resistance_genes.csv
    """
}

process callVirulenceARGENT {
    tag { project }

    publishDir "${outDir}/${project}/${RunID}/analysis", pattern: "*_virulence_genes.csv", mode: 'copy'

    container "file:///${params.simgdir}/abricate.simg"

    input:
    set project, file(assembly) from virAssemblyARGENT.groupTuple(by: 0)

    output:
    file "*_virulence_genes.csv"

    script:
    """
    abricate --db ${params.virulencedb} --csv *.fna > ${RunID}_virulence_genes.csv
    """
}


// Setup DIGCD trimmed reads channel

TrimmedReadsDIGCD.filter { it[1] ==~ 'digcd' }.into{ DIGCDTrimmedReadsShovill; DIGCDTrimmedReadsSnapperDB }

DIGCDMashRef = Channel.fromPath( "${params.digcdmashref}" )

DIGCDSnapperDBConfDir = Channel.fromPath( "${params.digcdsnapperdbconfdir}", type: 'dir', maxDepth: 1)

DIGCDSnapperDBRefDir = Channel.fromPath( "${params.digcdsnapperdbrefdir}", type: 'dir', maxDepth: 1)

float maxmashdist = Float.parseFloat(params.digcdmaxmashdist)

process assembleShovillDIGCD {
    tag { dataset_id }

    errorStrategy 'ignore'

    container "file:///${params.simgdir}/shovill.simg"

    cpus 4

    input:
    set dataset_id, project, file(forward), file(reverse) from DIGCDTrimmedReadsShovill

    output:
    set dataset_id, project, file("${dataset_id}.fasta") optional true into prokkaAssemblyDIGCD
    set project, file("${dataset_id}.fasta") optional true into quastDIGCD
    set dataset_id, project, file("${dataset_id}.fasta") into referenceSelection

    script:
    """
    shovill --cpus ${task.cpus} --R1 ${forward} --R2 ${reverse} --minlen 500 --outdir shovill
    mv shovill/contigs.fa ${dataset_id}.fasta
    """
}

process quastDIGCD {
    tag { project }

    publishDir "${outDir}/${project}/${RunID}/analysis", pattern: "*_assembly_summary.csv", mode: 'copy'

    container "file:///${params.simgdir}/quast.simg"

    cpus 4

    input:
    set project, file(assembly) from quastDIGCD.groupTuple(by: 0)

    output:
    file("*_assembly_summary.csv")

    script:
    """
    quast.py -t ${task.cpus} -o . --contig-thresholds 0 --no-html --no-plots *.fasta
    sed 's/\t/,/g' transposed_report.tsv > ${RunID}_assembly_summary.csv
    """
}

process annotateProkkaDIGCD {
    tag { dataset_id }

    publishDir "${outDir}/${project}/${RunID}/analysis/annotation/", pattern: "${dataset_id}/${dataset_id}.*", mode: 'copy'

    container "file:///${params.simgdir}/prokka.simg"

    cpus 4

    input:
    set dataset_id, project, file(assembly) from prokkaAssemblyDIGCD.filter{ it[2].size()>1000 }

    output:
    file("${dataset_id}/${dataset_id}.*")

    script:
    locusTag = dataset_id.tokenize("_")[0].tokenize("-")[1]
    prefix = project.toUpperCase()
    """
    prokka --outdir ${dataset_id} --locustag ${prefix}_${locusTag} --prefix ${dataset_id} --centre PHW --cpus ${task.cpus} --compliant ${assembly}
    cp ${dataset_id}/${dataset_id}.fna ${dataset_id}.fna
    """
}

process mashDistanceToRef {
    tag { dataset_id }

    publishDir "${outDir}/${project}/${RunID}/analysis/reference_selection/", pattern: "${dataset_id}_distance.csv", mode: 'copy'

    container "file:///${params.simgdir}/mash.simg"

    cpus 4

    input:
    set dataset_id, project, file(assembly), file(ref) from referenceSelection.combine(DIGCDMashRef).filter{ it[2].size()>1000 }

    output:
    set dataset_id, stdout into snapperClosestRef
    file("${dataset_id}_distance.csv")
    set project, file("${dataset_id}_closest.csv") into mashDistanceSummary
    set project, dataset_id, file("${dataset_id}_closest.csv") into createDistanceFilter

    script:
    """
    echo "#Reference,Query,Distance,P-value,Common kmers" > header.csv
    mash dist ${ref} ${assembly} | sort -k3 | sed 's/\t/,/g' | tee distance.csv | head -n1 | tee ${dataset_id}_closest.csv | cut -d "," -f1 | tr -d '\n'
    cat header.csv distance.csv > "${dataset_id}_distance.csv"
    """
}

process createDistanceFilter {
    tag { dataset_id }

    executor 'local'

    input:
    set project, dataset_id, file(tophit) from createDistanceFilter

    output:
    set dataset_id, stdout into distanceFilter

    script:
    """
    cut -d "," -f3 ${tophit} | tr -d '\n'
    """
}

process mashDistanceSummary {
    tag { RunID }

    publishDir "${outDir}/${project}/${RunID}/analysis/", pattern: "${RunID}_distance_summary.csv", mode: 'copy'

    executor 'local'

    input:
    set project, file(mashclosest) from mashDistanceSummary.groupTuple(by: 0)

    output:
    file("${RunID}_distance_summary.csv")


    script:
    """
    echo "#Reference,Query,Distance,P-value,Common kmers" > ${RunID}_distance_summary.csv
    cat *_closest.csv | sort -k3 >> ${RunID}_distance_summary.csv
    """
}

snapperDBInputUnfiltered = snapperClosestRef.combine(distanceFilter, by:0).combine(DIGCDTrimmedReadsSnapperDB, by: 0).combine(DIGCDSnapperDBConfDir).combine(DIGCDSnapperDBRefDir)

snapperDBrenameReads = Channel.create()

snapperDBInputfilterFail = Channel.create()

snapperDBInputUnfiltered.choice( snapperDBrenameReads , snapperDBInputfilterFail ){ float dist = Float.parseFloat(it[2])
                                                                                          dist <= maxmashdist ? 0 : 1 }


process renameReadsForSnapperDB {
    tag { dataset_id }

    input:
    set dataset_id, ref, distance, project, file(forward), file(reverse), file(configdir), file(refdir) from snapperDBrenameReads

    output:
    set dataset_id, ref, distance, project, file("${snapperDBname}.R1.fq.gz"), file("${snapperDBname}.R2.fq.gz"), file("${configdir}"), file("${refdir}") into snapperDBInputfilterPass

    script:
    snapperDBname = dataset_id.replace("DIGCD-", "").replaceAll(/_S\d+_L001$/, "")
    """
    mv ${forward} ${snapperDBname}.R1.fq.gz
    mv ${reverse} ${snapperDBname}.R2.fq.gz
    """
}


process snapperDBfastqToVCF {
    tag { dataset_id }

    publishDir "${outDir}/${project}/${RunID}/analysis/variant_calling/", pattern: "*vcf*", mode: 'copy'

    container "file:///${params.simgdir}/snapperdb.simg"

    cpus 4

    input:
    set dataset_id, ref, distance, project, file(forward), file(reverse), file(configdir), file(refdir) from snapperDBInputfilterPass

    output:
    file("*vcf*")
    set dataset_id, ref, project, file("*.filtered.vcf"), file("*.vcf.idx") into snapperDBVCFtoDB

    script:
    confname = ref.take(ref.lastIndexOf('.'))
    """
    export GASTROSNAPPER_CONFPATH=${configdir}
    export GASTROSNAPPER_REFPATH=${refdir}
    run_snapperdb.py fastq_to_vcf -c ${confname}.txt ${forward} ${reverse}
    mv snpdb/*.vcf* .
    """
}

/*
process snapperDBVCFtoDB {
   tag { dataset_id }

   maxForks 1

   cpus 1

   input:
   set dataset_id, ref, project, file(vcf), file(index), file(configdir), file(refdir) from snapperDBVCFtoDB

   output:
   val ref into snapperDBupdateDistanceMatrix
   val dataset_id into snapperDBupdateDistanceMatrixDummy

   script:
   confname = ref.take(ref.lastIndexOf('.'))
   """
   export GASTROSNAPPER_CONFPATH=${configdir}
   export GASTROSNAPPER_REFPATH=${refdir}
   run_snapperdb.py vcf_to_db -c ${confname}.txt ${vcf}
   """
}


process snapperDBupdateDistanceMatrix {
   tag { confname }

   cpus 1

   input:
   set ref, file(configdir) from snapperDBupdateDistanceMatrix.combine(confDirUpdateDistanceMatrix)
                                                              .combine(snapperDBupdateDistanceMatrixDummy.collect())
                                                              .unique()
                                                              .map{ [ it[0] , it[1] ] }

   output:
   val ref

   script:
   confname = ref.take(ref.lastIndexOf('.'))
   """
   export GASTROSNAPPER_CONFPATH=${configdir}
   run_snapperdb.py update_distance_matrix -c ${confname}.txt
   """
}
*/



workflow.onComplete {

    def msg = """\
        Pipeline execution summary
        ---------------------------
        MiSeq run name  : ${RunID}
        Workflow name   : ${workflow.scriptName}
        Workflow version: ${workflow.scriptId}
        Completed at    : ${workflow.complete}
        Duration        : ${workflow.duration}
        Success         : ${workflow.success}
        workDir         : ${workflow.workDir}
        exit status     : ${workflow.exitStatus}
        """
        .stripIndent()

if(workflow.success == true) {
    sendMail(to: 'matthew.bull@wales.nhs.uk', subject: "PenGU sequencing pipeline complete: SUCCESS", body: msg)
}
else {
    sendMail(to: 'matthew.bull@wales.nhs.uk', subject: "PenGU sequencing pipeline complete: FAILURE", body: msg)
}
}
	
