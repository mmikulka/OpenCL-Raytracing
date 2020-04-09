__kernel void NumericalReduction(__global int* data, __local int* localData, __global int *outData)
{
    size_t globalId = get_global_id(0);
    size_t localSize = get_local_size(0);
    size_t localId = get_local_id(0);
    
    localData[localId] = data[globalId];
    barrier(CLK_LOCAL_MEM_FENCE);
    
    for (int i = localSize >> 1; i > 0 ; i >>= 1)
    {
        if (localId < i)
            localData[localId] += localData[localId + i];
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    outData[get_group_id(0)] = 20;
}

