import socket
import sys

def main():
    HOST, PORT = "192.168.1.35", 8000
    data = " ".join(sys.argv[1:])

    # Create a socket (SOCK_STREAM means a TCP socket)
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

    try:
        # Connect to server and send data
        sock.connect((HOST, PORT))
        sock.sendall(data.encode())

        # Receive data from the server and shut down
        received = sock.recv(1024)
    finally:
        sock.close()

    print("Sent:     {}".format(data))
    print("Received: {}".format(received))

if __name__ == "__main__":
    main()