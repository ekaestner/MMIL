function x = empty2nan(x)

x(isempty(x)) = NaN;
return
end