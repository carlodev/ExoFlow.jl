using Interpolations, XLSX, DataFrames

function dns_velocity()
uτ = 5.7455e-2
ν = 1.4545e-4
U_DNS = DataFrame(XLSX.readtable("../src/Channel/Benchmark/Ret_395/U.xlsx", "Sheet1"))

y_physical = reverse(1 .- U_DNS.yp*ν/uτ)
u_physical = reverse(U_DNS.Umean*uτ)

u_DNS = Interpolations.LinearInterpolation(y_physical, u_physical, extrapolation_bc=Line())


function u_DNS_r(x::VectorValue)
   D = length(x)
   ux = u_DNS(abs(x[2]))
   if D == 2
      return VectorValue(ux, 0.0)
   elseif D == 3
      return VectorValue(ux, 0.0, 0.0)
   end
end


return u_DNS_r

end
