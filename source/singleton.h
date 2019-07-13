// Created by Ruiqi Zhang on 2019-07-09

#ifndef SINGLETON_H
#define SINGLETON_H

template<class T>
class Singleton
{
	public:
		static T& instance()
		{
			static T _instance;
			return _instance;
		}

	public:
		Singleton(const Singleton&) = delete;
		Singleton& operator=(const Singleton&) = delete;

	protected:
		Singleton() = default;
		virtual ~Singleton() = default;
};

#endif // SINGLETON_H
