#!/bin/bash

set -x

export SLURM_EXPORT_ENV='ALL'
#xport OMPI_MCA_pml=ucx

#xport DAPL_ACK_RETRY='7'
#xport DAPL_ACK_TIMER='20'
#xport DAPL_IB_SL='0'
#xport DAPL_MAX_INLINE='64'
#xport DAPL_RNR_RETRY='10'
#xport DAPL_RNR_TIMER='12'
#xport DAPL_UCM_CQ_SIZE='8192'
#xport DAPL_UCM_QP_SIZE='8192'
#xport DAPL_UCM_REP_TIME='10000'
#xport DAPL_UCM_RETRY='7'
#xport DAPL_UCM_RTU_TIME='8000'
#xport DAPL_UCM_TX_BURST='100'
#xport PMI_TIME='10'

exec $*

