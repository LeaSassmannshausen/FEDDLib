set(Checkpointing_HEADERS
    Checkpointing/CheckpointFiles.hpp
    Checkpointing/CheckpointMetadata.hpp
    Checkpointing/CheckpointMeshFingerprint.hpp
    Checkpointing/FSIOutletState.hpp
    Checkpointing/RecoveryCheckpoint.hpp
CACHE INTERNAL "")

set(LOCAL_HEADERS)
foreach(INC ${Checkpointing_HEADERS})
    string(REPLACE "Checkpointing/" "" LOCAL_INC ${INC})
    append_set(LOCAL_HEADERS ${LOCAL_INC})
endforeach()
