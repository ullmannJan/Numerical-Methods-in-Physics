using Plots


function f(x)
	return exp(sqrt(5)*x) - 13.5*cos(0.1*x)+25*x^4
end

function root_bisection(func, a, b, eps=1e-10)
	@assert a < b "a !< b"
	while b-a > eps
		mid = (a+b)/2
		if func(mid)*func(a) > 0
			a = mid
		else
			b = mid
		end
	end
    return a
end

function root_lin_interpol(func, a, b, eps=1e-14)
	@assert a < b "a !< b"
	new = 1
	count = 0
	while func(new) > eps || count > 1e4
		count++
		new = b-func(b)*(b-a)/(func(b)-func(a))
		if func(new)*func(a) > 0
			a = new
		else
			b = new
		end
	end
    return new
end


x = 1:10; y = rand(10); # These are the plotting data
plot(x,y)
# savefig("test.pdf")