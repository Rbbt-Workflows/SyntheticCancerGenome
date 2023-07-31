module SyntheticCancerGenome
  def self.choose(values)
    num = Misc.digest(values.inspect).hex[-1]
    i = num % values.length
    values[i]
  end

end
