import os

from parsl.config import Config
from parsl.executors import HighThroughputExecutor
from parsl.addresses import address_by_hostname
from parsl.launchers import SingleNodeLauncher
from parsl.providers import TorqueProvider


config = Config(
    executors=[
        HighThroughputExecutor(
            cores_per_worker=6,
            worker_debug=False,
            address=address_by_hostname(),
            provider=TorqueProvider(
                launcher=SingleNodeLauncher(),
                worker_init="export PYTHONPATH='{}:{{PYTHONPATH}}'".format(os.getcwd()),
                init_blocks=1,
                max_blocks=30,
                min_blocks=1,
                nodes_per_block=1,
                walltime='00:320:00'
            ),
        )
    ],
    retries=3
)
