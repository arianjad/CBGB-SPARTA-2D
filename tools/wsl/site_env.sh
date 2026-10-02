# Site-specific environment for the batch scripts; source it, then call
# site_mpi (SPARTA runs) or site_python (post-processing).
#   ORCD Engaging: gcc/12.2.0 + openmpi/4.1.4 modules; python 3.10 from the
#                  deprecated-modules tree, then the repo's .venv.
#   MIT subMIT:    `module load mpi` (central OpenMPI); the repo's .venv is
#                  built from a conda python 3.10, so no module is needed.
# Detection: subMIT mounts /ceph/submit; anything else is treated as ORCD.
if [ -d /ceph/submit ]; then SITE=submit; else SITE=orcd; fi
site_mpi() {
  if [ "$SITE" = submit ]; then module load mpi
  else module purge; module load gcc/12.2.0 openmpi/4.1.4; fi
}
site_python() {
  if [ "$SITE" = orcd ]; then module purge; module load deprecated-modules gcc/12.2.0-x86_64 python/3.10.8-x86_64; fi
}
