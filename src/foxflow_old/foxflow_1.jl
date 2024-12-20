using ZMQ
using JSON


"""
Dict{String, Any}(
    "shell_port" => 36203,
    "ip" => "127.0.0.1",
    "control_port" => 33947,
    "stdin_port" => 60695,
    "kernel_name" => "foxflow-0.1",
    "transport" => "tcp",
    "signature_scheme" => "hmac-sha256",
    "hb_port" => 47063,
    "jupyter_session" => "/home/rany/Work/playground/foxflow/Untitled3.ipynb",
    "key" => "9a7bf19e-740eafecde64e5fbd8c83587",
    "iopub_port" => 45999
)

"""

function load_connection_file(file_path)
    return JSON.parse(read(file_path, String))
end

function handle_shell_requests(shell_socket, iopub_socket)


    result_message = Dict(
    "header" => Dict(
        "msg_type" => "execute_result",
        "msg_id" => "unique_message_id",  # Generate a unique ID
        "username" => "username",          # Use a valid username
        "session" => "session_id",         # Use the correct session ID
        "version" => "5.0"
    ),
    "content" => Dict(
        "execution_count" => 1,  # Increment this for each execution
        "data" => Dict(
            "text/plain" => "Hello, World!",  # Your result
            "text/html" => "<div>Hello, World!</div>"
        ),
        "metadata" => Dict()
    )
)
    while true
        message = String(ZMQ.recv(shell_socket))
        println("Received shell request: $message")

        # Simulate some work
        sleep(1)

        # Send reply back to the client
        # ZMQ.send(shell_socket, "World shell")
        ZMQ.send(iopub_socket, JSON.json(handle_shell_requests))
    end
end

function handle_iopub_requests(iopub_socket)
    while true
        message = String(ZMQ.recv(iopub_socket))
        println("Received iopub message: $message")

        # Optionally process iopub messages
        ZMQ.send(iopub_socket, "World pub")
    end
end

function handle_message(shell_socket)

    try
        sleep(1)

        msg = recv(shell_socket)
        msg = JSON.parse(String(msg))

        msg_type = msg["header"]["msg_type"]
        
        if msg_type == "execute_request"
            code = msg["content"]["code"]
            result = eval(Meta.parse(code))
            send_execute_result(shell_socket, msg, result)
        end
    catch e
        rethrow(e)
    end
end

function send_execute_result(socket, msg, result)
    reply = Dict(
        "header" => Dict(
            "msg_type" => "execute_result",
            "msg_id" => msg["header"]["msg_id"],
            "username" => msg["header"]["username"],
            "session" => msg["header"]["session"],
            "version" => "5.0"
        ),
        "content" => Dict(
            "execution_count" => msg["content"]["execution_count"],
            "data" => Dict("text/plain" => string(result)),
            "metadata" => Dict()
        )
    )

    send(socket, JSON.json(reply))
end

function main()

    connection_file = ARGS[1]

    conn_info = load_connection_file(connection_file)

    context = ZMQ.Context()
    
    control_socket = ZMQ.Socket(context, ZMQ.REQ)
    stdin_socket = ZMQ.Socket(context, ZMQ.REQ)
    shell_socket = ZMQ.Socket(context, ZMQ.REQ)
    iopub_socket = ZMQ.Socket(context, ZMQ.PUB)
    hb_socket = ZMQ.Socket(context, ZMQ.PUB)

    println("conn_info: ", conn_info)

    # For handling control messages, such as interrupts.
    control_port = string(conn_info["control_port"])
    # For receiving input requests from the user.
    stdin_port = string(conn_info["stdin_port"])
    # For sending execution requests and receiving replies.
    shell_port = string(conn_info["shell_port"])
    # For sending output results, display data, and errors.
    iopub_port = string(conn_info["iopub_port"])
    # Heartbeat
    hb_port = string(conn_info["hb_port"])


    control_tcp = "tcp://*:$control_port"
    stdin_tcp = "tcp://*:$stdin_port"
    shell_tcp = "tcp://*:$shell_port"
    iopub_tcp = "tcp://*:$iopub_port"
    hb_tcp = "tcp://*:$hb_port"
    shell_tcp = "tcp://*:$shell_port"

    println("shell_tcp: ", shell_tcp)

    # Bind sockets to specified ports
    bind(control_socket, control_tcp)
    bind(stdin_socket, stdin_tcp)
    bind(shell_socket, shell_tcp)
    bind(iopub_socket, iopub_tcp)
    bind(hb_socket, hb_tcp)

    # try
        # bind(shell_socket, shell_tcp)
    # catch e
        # println("already connected to socket")
        # exit(0)
    # end

    # @async handle_message()
    @async handle_shell_requests(shell_socket, iopub_socket)
    # @async handle_iopub_requests(iopub_socket)

    # Main loop
    while true

        # Wait for next request from client
        # message0 = String(ZMQ.recv(shell_socket))
        # message1 = String(ZMQ.recv(iopub_socket))


        # println("Received request m0: $message0")
        # println("Received request m1: $message1")

        # Do some 'work'
        sleep(1)

        # Send reply back to client
        # ZMQ.send(shell_socket, "World shell")
        # ZMQ.send(iopub_socket, "World pub")

    end

    # classy hit men always clean up when finish the job.
    ZMQ.close(socket)
    ZMQ.close(context)

end

main()
