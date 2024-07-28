#pragma once
#include <cstring> 
#include <iostream> 
#include <winsock2.h>
#include <ws2tcpip.h>
#include <vector>


class DataTransfer
{
private:
	int port;
	SOCKET serverSocket;
	SOCKET clientSocket;
	sockaddr_in serverAddress;
	WSADATA wsaData;

	bool initializeWinsock();
	bool createSocket();
	bool bindSocket();
	bool listenForConnections();
	bool acceptConnection();
	void closeSocket();

public:
	DataTransfer(int port);
	~DataTransfer();

	bool start();
	void sendData(const std::unique_ptr<std::vector<std::vector<double>>>& data);

};
