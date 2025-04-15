process FASTTREE {
    label 'process_medium'

    container "quay.io/biocontainers/fasttree:2.1.10--h14c3975_3"

    input:
    path(aln)

    output:
    path("unrooted.fasttree.newick"), emit: treeGTR

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    OMP_NUM_THREADS=${task.cpus} FastTree -nt \\
        -gtr -gamma -spr 4 -mlacc 2 -slownni \\
        -out unrooted.fasttree.newick \\
        ${aln}
    """

    stub:
    def args = task.ext.args ?: ''
    
    """
    """
}
