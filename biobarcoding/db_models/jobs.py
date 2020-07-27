

# ALGORITHMS

class AlgorithmType:
    """ """
    pass


class Algorithm:
    """ Description of a specific algorithm"""
    pass


class Workflow:
    pass


class WorkflowStep:
    pass


class ComputeResourceProcessorArchitecture:
    pass


class ComputeResourceOS:
    pass


class ComputeResource:
    pass


class WorkflowInComputeResource:
    """
    Specific adaptation of a Worflow to a ComputeResource
    """
    pass


class Job:
    """ Core entity in the ALgorithmic domain.
    Workflows are launched into ComputeResources with some parameters and under an Identity
    Then, jobs can be checked, cancelled and results obtained
    Also, listing of filtered/ordered jobs may be retrieved
    """
    pass


class ProcessorArchitecture:
    pass


class ComputingArchitecture:
    pass
