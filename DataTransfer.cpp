#include "DataTransfer.h"

bool DataTransfer::initializeWinsock()
{
	return WSAStartup(MAKEWORD(2, 2), &wsaData) == 0;
}

bool DataTransfer::createSocket()
{
	this->serverSocket = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);
	return this->serverSocket != INVALID_SOCKET;
}

bool DataTransfer::bindSocket()
{
	this->serverAddress.sin_family = AF_INET;
	this->serverAddress.sin_addr.s_addr = INADDR_ANY;
	this->serverAddress.sin_port = htons(this->port);

	return bind(this->serverSocket, (sockaddr*)&this->serverAddress, sizeof(this->serverAddress)) != SOCKET_ERROR;
}

bool DataTransfer::listenForConnections()
{
	return listen(this->serverSocket, SOMAXCONN) != SOCKET_ERROR;
}

bool DataTransfer::acceptConnection()
{
	int addressLength = sizeof(this->serverAddress);
	this->clientSocket = accept(this->serverSocket, (sockaddr*)&this->serverAddress, &addressLength);

	return this->clientSocket != INVALID_SOCKET;
}

void DataTransfer::closeSocket()
{
	if (this->clientSocket != INVALID_SOCKET)
	{
		closesocket(this->clientSocket);
	}
	if (this->serverSocket != INVALID_SOCKET)
	{
		closesocket(this->serverSocket);
	}
}

DataTransfer::DataTransfer(int port) : port(port), serverSocket(INVALID_SOCKET), clientSocket(INVALID_SOCKET)
{

}

DataTransfer::~DataTransfer()
{
	closeSocket();
	WSACleanup();
}

bool DataTransfer::start()
{
	if (!initializeWinsock()) {
		std::cerr << "WSAStartup failed: " << WSAGetLastError() << std::endl;
		return false;
	}
	if (!createSocket()) {
		std::cerr << "Socket creation failed: " << WSAGetLastError() << std::endl;
		return false;
	}
	if (!bindSocket()) {
		std::cerr << "Bind failed: " << WSAGetLastError() << std::endl;
		return false;
	}
	if (!listenForConnections()) {
		std::cerr << "Listen failed: " << WSAGetLastError() << std::endl;
		return false;
	}
	if (!acceptConnection()) {
		std::cerr << "Accept failed: " << WSAGetLastError() << std::endl;
		return false;
	}
	return true;
}

void DataTransfer::sendData(const std::unique_ptr<std::vector<std::vector<double>>>& data)
{
	for (const auto& vec : *data)
	{
		size_t size = vec.size();
		send(this->clientSocket, reinterpret_cast<const char*>(&size), sizeof(size), 0); // Send the size of the vector first
		send(this->clientSocket, reinterpret_cast<const char*>(vec.data()), size * sizeof(double), 0);
	}
}


	
