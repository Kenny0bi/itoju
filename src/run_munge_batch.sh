#!/bin/bash
# Resumable batch: stream, filter and munge every endpoint in the definition families.
cd "$(dirname "$0")/.."
ENDPOINTS="G6_EPLEPSY FE FE_STRICT FE_MODE GE GE_MODE G6_STATUSEPI \
SLEEP G6_SLEEPAPNO G6_SLEEPAPNO_INCLAVO F5_INSOMNIA KRA_PSY_SLEEP_NONORG_EXMORE G6_SLEEPDISOTH F5_SLEEP_NOS \
K11_CONSTIPATION K11_OTHFUNC K11_IBS K11_FUNCDYSP K11_REFLUX \
N14_NEUROMUSCDYSBLADD N14_OTHBLADD \
F5_ADHD KRA_PSY_HYPERKIN_EXMORE F5_MILDRET KRA_PSY_MENTALRET_EXMORE KRA_PSY_DEVWIDE_EXMORE KRA_PSY_AUTISM_EXMORE"
for e in $ENDPOINTS; do
  [ -f "data/munged/$e.done" ] && continue
  nice -n 10 python3 src/01_fetch_munge_finngen.py "$e" >> logs/munge_batch.log 2>&1 || echo "FAILED $e $(date)" >> logs/munge_batch.log
done
echo "BATCH COMPLETE $(date)" >> logs/munge_batch.log
touch logs/munge_batch.complete
