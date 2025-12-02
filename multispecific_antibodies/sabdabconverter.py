print("Paste your Thera-SAbDab data, then press Enter twice when done:")
print("-" * 60)

# Read until double newline
lines = []
empty_count = 0

while empty_count < 2:
    try:
        line = input()
        if line.strip() == "":
            empty_count += 1
        else:
            empty_count = 0
            lines.append(line)
    except EOFError:
        break

text = '\n'.join(lines)

# Parse - handle case where labels and values are concatenated
import re

therapeutic = ""
targets = []
current_target = None

# Split on known keywords
parts = re.split(r'(Therapeutic|Target \d+|Heavy Chain \d+|Light Chain \d+)', text)

for i in range(len(parts)):
    part = parts[i].strip()
    
    if part == 'Therapeutic':
        if i+1 < len(parts):
            therapeutic = re.split(r'Target \d+', parts[i+1])[0].strip()
    
    elif re.match(r'Target \d+', part):
        current_target = {'name': '', 'hc': '', 'lc': ''}
        targets.append(current_target)
        if i+1 < len(parts):
            # Extract target name (everything before "Heavy Chain")
            target_text = re.split(r'Heavy Chain \d+', parts[i+1])[0].strip()
            current_target['name'] = target_text.replace('/', '-').replace(' ', '-')
    
    elif re.match(r'Heavy Chain \d+', part) and current_target:
        if i+1 < len(parts):
            # Extract sequence (everything before next keyword or structure info)
            seq_text = parts[i+1].split('Light Chain')[0].split('100%')[0].strip()
            current_target['hc'] = seq_text
    
    elif re.match(r'Light Chain \d+', part) and current_target:
        if i+1 < len(parts):
            # Extract sequence (everything before structure info)
            seq_text = parts[i+1].split('100%')[0].strip()
            current_target['lc'] = seq_text

# Build FASTA output
fasta_lines = []
for idx, t in enumerate(targets, 1):
    if t['hc']:
        fasta_lines.append(f">{therapeutic}_{t['name']}_Heavy_Chain_{idx}")
        fasta_lines.append(t['hc'])
    if t['lc']:
        fasta_lines.append(f">{therapeutic}_{t['name']}_Light_Chain_{idx}")
        fasta_lines.append(t['lc'])

fasta_output = '\n'.join(fasta_lines)

# Print to terminal
print("\n" + "=" * 60)
print("FASTA FORMAT OUTPUT:")
print("=" * 60 + "\n")
print(fasta_output)

# Create Python file
if therapeutic:
    filename = f"{therapeutic.lower()}.py"
    with open(filename, 'w') as f:
        f.write(f'{therapeutic.lower()} = """\n')
        f.write(fasta_output)
        f.write('\n"""\n')
    print("\n" + "=" * 60)
    print(f"File created: {filename}")
    print("=" * 60)