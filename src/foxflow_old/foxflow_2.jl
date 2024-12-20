using ZMQ
using JSON
using Printf


function load_connection_file(file_path)
    return JSON.parse(read(file_path, String))
end

function handle_shell_requests(shell_socket, iopub_socket)

    message = String(ZMQ.recv(iopub_socket))

    println("handling shell")

    print(@sprintf("%s", message))

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

    while true
        message = String(ZMQ.recv(shell_socket))
        println("Received shell request: $message")

        # Simulate some work
        sleep(1)

        # Send reply back to the client
        # ZMQ.send(shell_socket, "World shell")

        ZMQ.send(iopub_socket, JSON.json(message))
    end
end

function handle_iopub_requests(iopub_socket)
    while true
        message = String(ZMQ.recv(iopub_socket))
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

        println("Received iopub message: $message")

        # Optionally process iopub messages
        # ZMQ.send(iopub_socket, "World pub")

	ZMQ.send(iopub_socket, JSON.json(result_message))
    end
end

function heartbeat(hb_socket)
    while true
        ZMQ.send(hb_socket, "heartbeat")
        sleep(5)  # Send a heartbeat every 5 seconds
    end
end



function main()
    connection_file = ARGS[1]
    conn_info = load_connection_file(connection_file)

    hb_port = string(conn_info["hb_port"])
    hb_tcp = "tcp://*:$hb_port"

    context = ZMQ.Context()

    shell_socket = ZMQ.Socket(context, ZMQ.REQ)
    iopub_socket = ZMQ.Socket(context, ZMQ.PUB)
    hb_socket = ZMQ.Socket(context, ZMQ.PUB)

    bind(shell_socket, "tcp://*:" * string(conn_info["shell_port"]))
    bind(iopub_socket, "tcp://*:" * string(conn_info["iopub_port"]))
    bind(hb_socket, hb_tcp)

    # Start handling requests in separate async tasks
    @async handle_shell_requests(shell_socket, iopub_socket)
    @async heartbeat(hb_socket)
    # @async handle_iopub_requests(iopub_socket)
    # Start heartbeat in a separate task
    # @async heartbeat(iopub_socket)

    # Keep the main task alive
    while true
        println("waiting")
        sleep(1)  # This can be a heartbeat or simply keep the kernel running
    end
end

main()
