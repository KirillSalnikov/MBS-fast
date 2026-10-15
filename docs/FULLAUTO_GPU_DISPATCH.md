# Moving a CUDA checkpoint between cards

For an existing result, keep its original command, output path and `--fullauto-gpus` list. Add `--fullauto-resume-compatible OLD_BINARY --fullauto-dispatch-gpus 0` to run its saved work on visible GPU 0 only. The compatible binary must be the supported original executable that created the checkpoint. The CUDA precision, compiler and backend profile must match.

`--fullauto-dispatch-gpus` accepts distinct visible indices (1–8 devices). It changes worker placement and concurrency only. It leaves the saved estimator, execution recipe, sampling, optical settings and cache identities unchanged. `fullauto_parameters.json` records the original `selected_gpus` and the actual `dispatch_gpus` and `active_worker_slots`; `worker_start` events record the actual card. Both process and persistent workers obey the override.

Use ordinary `--fullauto-gpus 0` for a new single-card calculation. Dispatch override requires compatible resume. Stop the former controller before starting its replacement. Never run two controllers on the same output directory. Changing cards does not accelerate each task by itself; independent tasks can overlap CPU preparation and GPU processing, and throughput must be measured.
