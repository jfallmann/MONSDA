workflow.onComplete {
    script:
    """
    echo 'Workflow finished, no error'
    """
}
