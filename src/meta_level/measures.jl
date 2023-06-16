function KGL01(S::T) where {T<:NamedTuple}
  return (S.shared+S.right+S.left)/((S.shared + (S.shared+S.right+S.left))/2.0)
end

function KGL02(S::T) where {T<:NamedTuple}
  return KGL01(S).-1.0
end

function KGL03(S::T) where {T<:NamedTuple}
  return (S.right + S.left)/2.0
end

function KGL04(S::T) where {T<:NamedTuple}
  return S.right + S.left
end

function KGL05(S::T) where {T<:NamedTuple}
  return ((S.shared+S.right+S.left)^2)/((S.shared+S.right+S.left)^2-2*S.right*S.left)-1.0
end

function KGL06(S::T) where {T<:NamedTuple}
  one_over_sum = 1/(2*S.shared+S.right+S.left)
  s1 = (S.shared+S.right)*log(S.shared+S.right)
  s2 = (S.shared+S.left)*log(S.shared+S.left)
  return log(2*S.shared+S.right+S.left)-(one_over_sum*2*S.shared*log(2))-(one_over_sum*(s1+s2))
end

function KGL07(S::T) where {T<:NamedTuple}
  return exp(KGL06(S))-1.0
end

function KGL08(S::T) where {T<:NamedTuple}
  return (S.right + S.left) / (2*S.shared + S.right + S.left)
end

function KGL09(S::T) where {T<:NamedTuple}
  return KGL08(S)
end

function KGL10(S::T) where {T<:NamedTuple}
  return S.shared / (S.shared + S.right + S.left)
end

function KGL11(S::T) where {T<:NamedTuple}
  return (2*S.shared) / (2*S.shared + S.right + S.left)
end

function KGL12(S::T) where {T<:NamedTuple}
  return (2*S.shared+S.right+S.left) * ( 1 - (S.shared/(S.shared+S.right+S.left)))
end

function KGL13(S::T) where {T<:NamedTuple}
  return min(S.right, S.left) / (max(S.right, S.left) + S.shared)
end

function KGL14(S::T) where {T<:NamedTuple}
  return 1 - (S.shared*(2*S.shared+S.right+S.left))/(2*(S.shared+S.right)*(S.shared+S.left))
end

function KGL15(S::T) where {T<:NamedTuple}
  return (S.right + S.left)/(S.shared + S.right + S.left)
end

function KGL16(S::T) where {T<:NamedTuple}
  return KGL15(S)
end

function KGL17(S::T) where {T<:NamedTuple}
  return min(S.right, S.left)/(S.shared+S.right+S.left)
end

function KGL18(S::T) where {T<:NamedTuple}
  return (S.right+S.left)/2
end

function KGL19(S::T) where {T<:NamedTuple}
  return ((S.right*S.left)+1)/((S.shared+S.right+S.left)^2-(S.shared+S.right+S.left)/2)
end

function KGL20(S::T) where {T<:NamedTuple}
  return 1 - (2*S.shared)/(2*S.shared.+S.right+S.left)
end

function KGL21(S::T) where {T<:NamedTuple}
  return S.shared/(S.shared+S.left)
end

function KGL22(S::T) where {T<:NamedTuple}
  return min(S.right, S.left)/(min(S.right, S.left) + S.shared)
end

function KGL23(S::T) where {T<:NamedTuple}
  return 2*abs(S.right-S.left)/(2*S.shared+S.right+S.left)
end

function KGL24(S::T) where {T<:NamedTuple}
  in_par = (2*S.shared+S.right+S.left)/(S.shared+S.right+S.left)
  return 1-log(in_par)/log(2)
end
