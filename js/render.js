const render = (thet,P0,h,t1,t2,i) => {
    const app = document.querySelector('#app');
    let class_name = `data-${i}`;
    const container = document.createElement('div');
    container.classList.add('app__box');
    app.insertAdjacentElement('beforeend', container);
    const block = `
    <h2 class="app__title">Расчеты №${i}</h2>
    <ul class="app__fix"></ul>
    <div class="app__table">
    <table>
    <thead>
    <tr>
    <th>t, с</th>
    <th>m, кг</th>
    <th>Vx, м/с</th>
    <th>Vy, м/с</th>
    <th>x, м</th>
    <th>y, м </th>
    <th>V, м/с </th>
    <th>r, м</th>
    <th>thet, рад</th>
    <th>THET , рад</th>
    <th>alfa, рад</th>
    <th>fi, рад</th>
    <th>THETc, рад</th>
    </tr>
    </thead>
    </table>
    <div class="app__table-body">
    <table>
    <tbody id="${class_name}">
    </tbody>
    </table>
    </div>
    </div>`
    container.insertAdjacentHTML('beforeend', block);
    const fixDiv = container.querySelector('.app__fix');
    
    fixDiv.innerHTML = `
    <li><b>P = </b>${P0} [Н]</li>
    <li><b>h = </b>${h} [м]</li>
    <li><b>t1 = </b>${t1} [c]</li>
    <li><b>t2 = </b>${t2} [c]</li>
    <li><b>thet_torch = </b>${thet[0]} [рад]</li>
    <li><b>thet_2 = </b>${thet[1]} [рад]</li>
    `
}
const html = (t, el, i) => {
    let class_name = `data-${i}`;
    const tr_iter = (t, el) => {
        return `
        <td class="app__ceil">${t.toFixed(5)}</td>
        <td class="app__ceil">${el[0].toFixed(5)}</td>
        <td class="app__ceil">${el[1].toFixed(5)}</td>
        <td class="app__ceil">${el[2].toFixed(5)}</td>
        <td class="app__ceil">${el[3].toFixed(5)}</td>
        <td class="app__ceil">${el[4].toFixed(5)}</td>
        <td class="app__ceil">${el[5].toFixed(5)}</td>
        <td class="app__ceil">${el[6].toFixed(5)}</td>
        <td class="app__ceil">${el[7].toFixed(5)}</td>
        <td class="app__ceil">${el[8].toFixed(5)}</td>
        <td class="app__ceil">${el[9].toFixed(5)}</td>
        <td class="app__ceil">${el[10].toFixed(5)}</td>
        <td class="app__ceil">${el[11].toFixed(5)}</td>
        `
    }
    
    const tbody = document.querySelector(`#${class_name}`);
    const tr = document.createElement('tr')
    tr.insertAdjacentHTML('beforeend', tr_iter(t, el))
    tbody.insertAdjacentElement('beforeend', tr)
}

const consoleLog = (t, el) => {
    console.log(`${t.toFixed(4)} ${el[0].toFixed(4)} ${el[1].toFixed(4)} ${el[2].toFixed(4)} ${el[3].toFixed(4)} ${el[4].toFixed(4)} ${el[5].toFixed(4)} ${el[6].toFixed(4)} ${el[7].toFixed(4)} ${el[8].toFixed(4)} ${el[9].toFixed(4)} ${el[10].toFixed(4)} ${el[11].toFixed(4)}`);
}

export {
    html,
    render
}