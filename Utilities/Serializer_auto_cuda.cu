#include "Serializer.h"
#include <stdio.h>

BEGIN_CUDA_COMPATIBLE()

CSerializer::CSerializer()
{
    readPos_ = 0;
}

void CSerializer::operator<<(int i)
{
    archive_.emplace_back(std::to_string(i));
}

void CSerializer::operator<<(size_t i)
{
    archive_.emplace_back(std::to_string(i));
}

void CSerializer::operator<<(char c)
{
    archive_.emplace_back(std::to_string(c));
}

void CSerializer::operator<<(double f)
{
    archive_.emplace_back(std::to_string(f));
}

void CSerializer::operator<<(float f)
{
    archive_.emplace_back(std::to_string(f));
}

void CSerializer::operator<<(bool b)
{
    archive_.emplace_back(std::to_string(b ? 1 : 0));
}

void CSerializer::operator<<(std::string s)
{
    archive_.emplace_back(s);
}

void CSerializer::operator>>(int& i)
{
    i = (int)std::atoi(archive_[readPos_++].data());
}

void CSerializer::operator>>(size_t& i)
{
    i = (size_t)std::atoi(archive_[readPos_++].data());
}

void CSerializer::operator>>(char& c)
{
    c = (char)std::atoi(archive_[readPos_++].data());
}

void CSerializer::operator>>(double& f)
{
    f = (double)std::atof(archive_[readPos_++].data());
}

void CSerializer::operator>>(float& f)
{
    f = (float)std::atof(archive_[readPos_++].data());
}

void CSerializer::operator>>(bool& b)
{
    b = (std::atoi(archive_[readPos_++].data()) == 1) ? true : false;
}

void CSerializer::operator>>(std::string& s)
{
    s = archive_[readPos_++];
}

void CSerializer::seekBegin()
{
    readPos_ = 0;
}

END_CUDA_COMPATIBLE()
