function main()
    a = 0
    for i in 1:3000000000
        a += i^2
        a /= 2
        a += 52
    end
    show(a)
end
@time main()
