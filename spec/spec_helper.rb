THIS_DIR = File.dirname __FILE__

require_relative File.join "..", "lib", "lsa.rb"

Dir[File.join("..", "lib", "lsa", "*")].each do |fname|
  require_relative File.join fname
end

RSpec.configure do |config|
  config.expect_with :rspec do |c|
    c.syntax = :expect
  end
end
