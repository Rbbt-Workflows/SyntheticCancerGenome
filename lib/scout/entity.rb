require 'scout/meta_extension'
module Entity
  extend MetaExtension

  extension_attr :entity_directory, :entity_type_directory, :properties, :workflow

  def job_options(task, options = {})
    options 
  end

  def job(task, options = {} )
    jobname = self
    job_options = job_options(task, jobname, options)
    self.job(job_options, jobname, job_options)
  end

  def define_property(property, options = {}, &block)
    @properties ||= {}
    @properties[property] = block
  end

  def exec_property(property, entity = nil, options = {})
    @properties[property].call entity, options
  end

  def property(property, entity = nil, options = {}, &block)
    return define_property(property, options, &block) if block_given?
    @cache ||= {}
    @@cache ||= {}
    return @cache[property] if @cache.include?(property)
    return @@cache[property][entity] if entity && @@cache[property] && @@cache[property].include?(entity)
    begin
      return @result = entity_directory[property] if entity_directory && entity_directory[property].exists?
      return @result = entity_type_directory[property][entity] if entity && entity_type_directory && entity_type_directory[property][entity].exists?
      return @result = entity_type_directory[property] if entity_type_directory && entity_type_directory[property].exists?
      return @result = exec_property(property, entity, options) if properties && properties.include?(property)
      return @result = workflow.job(property, entity, options) if workflow && workflow.tasks.include?(property)
      return @result = self.send(property) if self.respond_to?(property)
    ensure
      @cache[property] = @result
      @@cache[property] = @result
    end
  end

  def file(...)
    res = property(...)
    res
  end
end
