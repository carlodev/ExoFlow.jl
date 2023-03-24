
const PArrays=PartitionedArrays

function mytic!(t,comm)
    MPI.Barrier(comm)
    PArrays.tic!(t)
end