#ifndef AUTO_POINTER_INCLUDED
#define AUTO_POINTER_INCLUDED

template<class Type>
class auto_ptr 
{
public:
	typedef Type element_type;
	explicit auto_ptr(Type* ptr = 0) throw() : _own(ptr != 0), _pointer(ptr) 
	{
	}

	auto_ptr(const auto_ptr<Type>& _Y) throw() : _own(_Y._own), _pointer(_Y.release()) 
	{
	}

	auto_ptr<Type>& operator=(const auto_ptr<Type>& _Y) throw()
	{
		if (_pointer != _Y.get())
		{
			if (_own)
				delete _pointer;
			_own = _Y._own;
			_pointer = _Y.release(); 
		}
		else 
			if (_Y._own)
				_own = true;
		return *this; 
	}
	
	~auto_ptr()
	{
		if (_own)
			delete _pointer; 
	}

	Type& operator*() const throw()
	{
		return *get(); 
	}

	Type *operator->() const throw()
	{
		return get(); 
	}
	
	Type *get() const throw()
	{
		return _pointer; 
	}

	Type *release() const throw()
	{
		_own = false;
		return _pointer; 
	}
private:
	mutable bool _own;
	Type* _pointer;
};

