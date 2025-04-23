Misc.add_libdir if __FILE__ == $0

#require 'scout/sources/SyntheticCancerGenome'

Workflow.require_workflow "HTS"
Workflow.require_workflow "NEATGenReads"

module SyntheticCancerGenome
  extend Workflow

end

Workflow.main = SyntheticCancerGenome

#require 'synthetic_cancer_genome/tasks/somatic.rb'
#require 'synthetic_cancer_genome/tasks/minify.rb'

#require 'synthetic_cancer_genome/tasks/evolution.rb'
#

require 'synthetic_cancer_genome/tasks/clonality.rb'
require 'synthetic_cancer_genome/tasks/bundle.rb'
require 'synthetic_cancer_genome/tasks/analyze.rb'
require 'synthetic_cancer_genome/tasks/germline.rb'

#require 'synthetic_cancer_genome/tasks/reference.rb'
#require 'synthetic_cancer_genome/tasks/samples.rb'
#require 'synthetic_cancer_genome/tasks/reference.rb'
#require 'scout/knowledge_base/SyntheticCancerGenome'
#require 'scout/entity/SyntheticCancerGenome'
