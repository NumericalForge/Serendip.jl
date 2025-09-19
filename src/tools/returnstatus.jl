mutable struct ReturnStatus
    successful::Bool
    message::String
    function ReturnStatus(successful::Bool=true, message::String="")
        return new(successful, message)
    end
end

failed(rs::ReturnStatus)    = !rs.successful
succeeded(rs::ReturnStatus) =  rs.successful

failure(msg...)      = ReturnStatus(false, join(msg))
Base.success(msg...) = ReturnStatus(true, join(msg))
