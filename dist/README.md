# Linux CUDA builds for c0c2d5ecac57

Source commit `c0c2d5ecac57adc6ef58637cbf178203a1a646ef`. Choose the build for the listed host or matching CPU, CUDA architecture and runtime.

| Hosts | File | CUDA arch | Precision |
|---|---|---|---|
| epyc1, epyc2 | [mbs-fast-epyc1-sm70-double.gz](mbs-fast-epyc1-sm70-double.gz) | sm_70 | double |
| gpu1 | [mbs-fast-gpu1-sm86-float.gz](mbs-fast-gpu1-sm86-float.gz) | sm_86 | float |
| gpu2 | [mbs-fast-gpu2-sm86-float.gz](mbs-fast-gpu2-sm86-float.gz) | sm_86 | float |
| dima, ilya | [mbs-fast-dima-sm86-float.gz](mbs-fast-dima-sm86-float.gz) | sm_86 | float |
| v100 | [mbs-fast-v100-sm70-double.gz](mbs-fast-v100-sm70-double.gz) | sm_70 | double |

Run `sha256sum -c SHA256SUMS --ignore-missing` after downloading. Unpack the selected file with `gzip -dc FILE.gz > mbs-fast`, then run `chmod +x mbs-fast` and `./mbs-fast --version`.
