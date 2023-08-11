process WRITE_NETWORK_TO_DIR {
    publishDir params.publish_dir

    input:
    val (tool)
    tuple val (key), path (network)

    output:
    path ("networks/${key}-${tool}-network.txt")

    script:
    """
    mkdir -p "networks/"
    touch "networks/${key}_${tool}_network.txt"
    cat "$network" >> "networks/${key}-${tool}-network.txt"
    """
}