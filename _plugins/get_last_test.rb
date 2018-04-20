
module Jekyll
  module TestLookups
    #! we send the reports from the site in the template 
    #! ... because we are not sure how to get site here?
    def get_last_test(input,reports)
      matches = []
      for entry in reports
        if entry.key?('validations')
          for subtype in ['automatic','interface']
            if entry['validations'].key?(subtype)
              for item in entry['validations'][subtype]
                if item['tag']==input
                  matches << entry['when']
                end
              end
            end
          end
        end
      end
      if matches.any?
        report_code = matches.sort.last
        out = "<a class=\"bubble_light\" href=\"/reports/#"+
          report_code+"\">tested "+report_code+"</a>"
        return out
      else
        return "<span class=\"bubble_light\">pending test</span>"
      end
    end
  end
end

Liquid::Template.register_filter(Jekyll::TestLookups)
