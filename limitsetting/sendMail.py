# Import smtplib for the actual sending function
import smtplib

# Import the email modules we'll need
from email.mime.text import MIMEText

# Open a plain text file for reading.  For this example, assume that
# the text file contains only ASCII characters.
#fp = open(textfile, 'rb')
# Create a text/plain message
#msg = MIMEText(fp.read())
#fp.close()

msg = MIMEText('tests')

# me == the sender's email address
# you == the recipient's email address
msg['Subject'] = 'Test'
msg['From'] = 'james.keaveney@vub.ac.be'
msg['To'] = 'keaveneyenator@gmail.com'

# Send the message via our own SMTP server, but don't include the
# envelope header.
s = smtplib.SMTP('localhost')
s.sendmail(me, [you], msg.as_string())
s.quit()
