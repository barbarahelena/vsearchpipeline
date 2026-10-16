nextflow.enable.types = true

process GET_GIT_COMMIT {
    label 'process_single_low'
    label 'git'

    input:
    commit_hash: String

    output:
    commit_hash = file("git_commit_hash.txt")

    topic:
    tuple(task.process, 'git', eval('git --version | cut -d" " -f3')) >> 'versions'

    script:
    """
    echo "\nYou can obtain the exact version of the vsearchpipeline used for this run" >> git_commit_hash.txt
    echo "by cloning the repo and running the Git command:\n" >> git_commit_hash.txt
    echo "git checkout ${commit_hash}\n" >> git_commit_hash.txt
    """

    stub:
    """
    touch git_commit_hash.txt
    """
}
