using ZMQ
using JSON
using Printf


struct JupyterMessageHeader
    msg_id::String
    session::String
    username::String
    date::String
    msg_type::String
end

struct JupyterMessage
    message::JupyterMessageHeader
    parent_header::JupyterMessageHeader
end

function (msg::JupyterMessage)()
    result_message = Dict
end

function load_connection_file(file_path)
    return JSON.parse(read(file_path, String))
end

function handle_shell_requests(shell_socket, iopub_socket)
    while true
        message = String(ZMQ.recv(shell_socket))
        # println("Received shell request: $message")
        println("Received shell request: $message")
        flush(stdout)

        # Simulate some work
        sleep(1)

        # Prepare a response message
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
                    "text/plain" => "Hello, World! $message",  # Your result
                    "text/html" => "<div>Hello, World!</div>"
                ),
                "metadata" => Dict()
            )
        )

        # Send the response back to the shell socket
        ZMQ.send(shell_socket, JSON.json(result_message))
        
        # Also publish the result to the iopub socket if needed
        ZMQ.send(iopub_socket, JSON.json(result_message))
    end
end

function heartbeat(hb_socket)
    while true
        ZMQ.send(hb_socket, "heartbeat")
        sleep(5)  # Send a heartbeat every 5 seconds
    end
end

function send_kernel_info(iopub_socket)
    info_message = Dict(
        "header" => Dict(
            "msg_type" => "kernel_info_reply",
            "msg_id" => "unique_message_id",  # Generate a unique ID
            "username" => "username",
            "session" => "session_id",
            "version" => "5.0"
        ),
        "content" => Dict(
            "protocol_version" => "5.0",
            "implementation" => "FoxFlow",
            "implementation_version" => "1.0",
            "language" => "foxflow",
            "language_version" => "0.1",
            "banner" => "FoxFlow Kernel"
        )
    )

    ZMQ.send(iopub_socket, JSON.json(info_message))
end

function main()

    connection_file = ARGS[1]
    conn_info = load_connection_file(connection_file)

    println("connection_file: ", connection_file)

    context = ZMQ.Context()

    shell_socket = ZMQ.Socket(context, ZMQ.REP)  # Change to REP
    iopub_socket = ZMQ.Socket(context, ZMQ.PUB)
    hb_socket = ZMQ.Socket(context, ZMQ.PUB)

    bind(shell_socket, "tcp://*:" * string(conn_info["shell_port"]))
    bind(iopub_socket, "tcp://*:" * string(conn_info["iopub_port"]))
    bind(hb_socket, "tcp://*:" * string(conn_info["hb_port"]))

    # Start handling requests in separate tasks
    @async handle_shell_requests(shell_socket, iopub_socket)
    @async heartbeat(hb_socket)

    println("iopub_socket")

    send_kernel_info(iopub_socket)

    # Keep the main task alive
    while true
        sleep(1)  # Keep the kernel running
    end

end

main()
