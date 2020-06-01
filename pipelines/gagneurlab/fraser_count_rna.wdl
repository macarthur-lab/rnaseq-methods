version 1.0

workflow FraserCountRNAWorkflow {

    input {
        File input_bam
        File input_bai
    }

    call FraserCountRNA {
        input:
            input_bam=input_bam,
            input_bai=input_bai
    }

    output {
        File output_cache_tar_gz = FraserCountRNA.output_cache_tar_gz
    }
}

task FraserCountRNA {

    input {
        File input_bam
        File input_bai

        String sample_id = sub(basename(input_bam), "\\.bam$|\\.cram$", "")

        Int disk_size = ceil(size(input_bam, "GB") + 15)
    }

    command {
        echo --------------; echo "Start - time: $(date)"; set -euxo pipefail; ls -lhtr; lscpu | grep 'CPU\|Model'; free -h; df -kh; uptime; find /cromwell*/ -type f | xargs ls -lhSr; echo --------------

        if [ ! -f ~{input_bam}.bai ]; then
            mv ~{input_bai} ~{input_bam}.bai
        fi

        time Rscript --vanilla /countRNA.R --num-threads 4 ~{sample_id} ~{input_bam}

        tar czf fraser_count_rna_~{sample_id}.tar.gz cache

        echo --------------; ls -lhtr; lscpu | grep 'CPU\|Model'; free -h; df -kh; uptime; set +xe; echo "Done - time: $(date)"; echo --------------
    }

    output {
        File output_cache_tar_gz = "fraser_count_rna_~{sample_id}.tar.gz"
    }

    runtime {
        docker: "weisburd/gagneurlab@sha256:dc3bdc7873f025c62fa8e1c1d4ae5e3670d4307baf8b8f24baae4e1bb39e3a80"
        cpu: 1
        preemptible: 1
        disks: "local-disk ${disk_size} HDD"
        #memory: "15 GB"
        #disks: "local-disk ${disk_size} LOCAL"
        #bootDiskSizeGb: 20
        #zones: "us-east1-b us-east1-c us-east1-d"
    }
}
