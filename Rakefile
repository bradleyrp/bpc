require 'html-proofer'

task :test do
  sh "bundle exec jekyll build"
  options = { :assume_extension => true, :empty_alt_ignore => true }
  HTMLProofer.check_directory("./_site", options).run
end
